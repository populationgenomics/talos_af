"""
whole new approach for testing MOI, very simple
"""

import abc
from collections import defaultdict
from dataclasses import dataclass
from typing import ClassVar

from mendelbrot.pedigree_parser import PedigreeParser

from talos_af import config, models
from talos_af import utils as utils_af
from talos_af.models import ReportableVariant

HEMI_CHROMS = {'chrX', 'chrY'}


@dataclass
class GlobalFilter:
    """
    A Filter class, used to apply to any non-ClinVar Pathogenic Variants
    This pulls in a number of thresholds from the config file, and contains a 'too_common' method
    A variant run through this method will return true based on the thresholds if it's 'too_common'
    """

    # a lookup of the attribute name vs. the corresponding configurable filter to be used on small variants
    small_dict: ClassVar[dict[str, float | int]] = {
        'gnomad_af': config.config_retrieve('gnomad_max_af'),
        'gnomad_homalt': config.config_retrieve('gnomad_max_homozygotes'),
    }

    # only to be applied on chrX/Y
    small_gnomad_hemi: ClassVar[int] = config.config_retrieve('gnomad_max_hemizygotes')

    def too_common(self, variant: models.VariantAf) -> bool:
        """Check if a variant is too common in the population."""
        for key, threshold in self.small_dict.items():
            if key in variant.info and variant.info[key] > threshold:
                return True
        # on sex chroms, apply hemi-count filter
        if variant.coordinates.chrom in HEMI_CHROMS:
            return variant.info.get('gnomad_ac_xy', 0) > self.small_gnomad_hemi
        return False


@dataclass
class DominantFilter:
    """
    Similar to the GlobalFilter, but with stricter thresholds
    This is designed to run on variants being considered for Dominant inheritance
    """

    # a lookup of the attribute name vs. the corresponding configurable filter
    small_dict: ClassVar[dict[str, float | int]] = {
        'gnomad_af': config.config_retrieve('dominant_gnomad_max_af'),
        'gnomad_ac': config.config_retrieve('dominant_gnomad_max_ac'),
        'gnomad_homalt': config.config_retrieve('dominant_gnomad_max_homozygotes'),
    }

    def too_common(self, variant: models.VariantAf) -> bool:
        """Check if a variant is too common in the population."""
        # check against each small-variant filter
        return any(key in variant.info and variant.info[key] > threshold for key, threshold in self.small_dict.items())


@dataclass
class ClinVarFilter:
    """This will apply more lenient filters to ClinVar Pathogenic variants."""

    # a lookup of the attribute name vs. the corresponding configurable filter
    small_dict: ClassVar[dict[str, float]] = {
        'gnomad_af': config.config_retrieve('clinvar_gnomad_max_af'),
    }

    def too_common(self, variant: models.VariantAf) -> bool:
        """Check if a variant is too common in the population."""
        return any(key in variant.info and variant.info[key] > threshold for key, threshold in self.small_dict.items())


@dataclass
class ClinVarDominantFilter:
    """This will apply more lenient filters to ClinVar Pathogenic variants, Designed to run on Dominant variants."""

    # a lookup of the attribute name vs. the corresponding configurable filter
    small_dict: ClassVar[dict[str, float]] = {
        'gnomad_af': config.config_retrieve('clinvar_dominant_gnomad_max_af'),
        'af': config.config_retrieve('clinvar_dominant_callset_max_af'),
    }

    def too_common(self, variant: models.VariantAf) -> bool:
        """Check if a variant is too common in the population."""
        return any(key in variant.info and variant.info[key] > threshold for key, threshold in self.small_dict.items())


class BaseMoi(abc.ABC):
    """
    Definition of the MOI test base class
    """

    def __init__(self, pedigree: PedigreeParser, applied_moi: str):
        if applied_moi is None:
            raise ValueError('An applied MOI needs to reach the Base Class')
        self.applied_moi = applied_moi
        self.pedigree = pedigree
        self.minimum_depth = config.config_retrieve('minimum_depth', 10)
        self.global_filter: GlobalFilter | DominantFilter = GlobalFilter()
        self.clinvar_filter: ClinVarFilter | ClinVarDominantFilter = ClinVarFilter()

    @abc.abstractmethod
    def run(
        self,
        principal: models.VariantAf,
        comp_het: utils_af.CompHetDict,
    ) -> dict[str, list[models.ReportableVariant]]:
        """
        run all applicable inheritance patterns and finds good fits
        """

    def variant_too_common(self, variant: models.VariantAf) -> bool:
        """
        Check if a variant is too common in the population or callset

        Args:
            variant (VariantAf): the variant to check

        Returns:
            bool: True if the variant is too common
        """
        if variant.info.get('clinvar_plp'):
            return self.clinvar_filter.too_common(variant=variant)
        return self.global_filter.too_common(variant)

    @staticmethod
    def comp_het_out_of_phase(
        sample_id: str,
        variant_1: models.VariantAf,
        variant_2: models.VariantAf,
    ) -> bool:
        """
        uses phase and genotypes to check if a comp-het is valid for a sample

        Args:
            sample_id (str): sample ID to check for
            variant_1 (models.VariantAf): first variant of comp-het pair
            variant_2 (models.VariantAf): second variant of comp-het pair

        Returns:
            bool: True if these two variants may form a comp-het
        """

        # check for sample id in the set intersection
        if sample_id not in (variant_1.het_samples & variant_2.het_samples):
            return False

        # check the two variants are out of phase
        # no phase data, no possible phase match
        if (sample_id not in variant_1.phased) or (sample_id not in variant_2.phased):
            return True

        for phase_set, genotype in variant_1.phased[sample_id].items():
            # this phase set doesn't exist in both variants - can't be phased
            if phase_set not in (v2p := variant_2.phased[sample_id]):
                return True
            return genotype != v2p

        return True


class AD(BaseMoi):
    def __init__(self, pedigree: PedigreeParser, applied_moi: str = 'Autosomal Dominant'):
        """
        Simplest: AD MOI
        """

        super().__init__(pedigree=pedigree, applied_moi=applied_moi)
        self.global_filter = DominantFilter()
        self.clinvar_filter = ClinVarDominantFilter()

    def run(
        self,
        principal: models.VariantAf,
        comp_het: utils_af.CompHetDict,  # noqa: ARG002
    ) -> dict[str, list[models.ReportableVariant]]:
        """Simplest MOI, exclusions based on HOM count and AF."""

        classifications: dict[str, list[models.ReportableVariant]] = defaultdict(list)

        if self.variant_too_common(principal):
            return classifications

        # autosomal dominant doesn't require support, but consider het and hom
        for sample_id in principal.het_samples.union(principal.hom_samples):
            classifications[sample_id].append(ReportableVariant(var_id=principal.coordinates.string_format))

        return classifications


class AR(BaseMoi):
    def __init__(self, pedigree: PedigreeParser, applied_moi: str = 'Autosomal Recessive'):
        """Autosomal Recessive MOI - check for homs and compound-hets."""

        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: models.VariantAf,
        comp_het: utils_af.CompHetDict,
    ) -> dict[str, list[models.ReportableVariant]]:
        """Simplest MOI, exclusions based on HOM count and AF."""

        classifications: dict[str, list[models.ReportableVariant]] = defaultdict(list)

        if self.variant_too_common(principal):
            return classifications

        # autosomal dominant doesn't require support, but consider het and hom
        for sample_id in principal.hom_samples:
            classifications[sample_id].append(ReportableVariant(var_id=principal.coordinates.string_format))

        for sample_id in principal.het_samples:
            support_vars = []
            for partner in comp_het[sample_id].get(principal.coordinates.string_format, []):
                if self.comp_het_out_of_phase(sample_id, principal, partner):
                    support_vars.append(partner)

            if support_vars:
                classifications[sample_id].append(
                    ReportableVariant(
                        var_id=principal.coordinates.string_format,
                        support_vars={partner.coordinates.string_format for partner in support_vars},
                    )
                )

        return classifications


class XL(BaseMoi):
    def __init__(self, pedigree: PedigreeParser, applied_moi: str = 'X-Linked'):
        """X-Linked - Heterozygous and homozygous variants in females, as well as hemizygous in males are reportable."""
        super().__init__(pedigree=pedigree, applied_moi=applied_moi)

    def run(
        self,
        principal: models.VariantAf,
        comp_het: utils_af.CompHetDict,  # noqa: ARG002
    ) -> dict[str, list[models.ReportableVariant]]:
        """Simplest MOI, exclusions based on HOM count and AF."""

        classifications: dict[str, list[models.ReportableVariant]] = defaultdict(list)

        if self.variant_too_common(principal):
            return classifications

        # autosomal dominant doesn't require support, but consider het and hom
        for sample_id in principal.hom_samples | principal.het_samples:
            classifications[sample_id].append(
                ReportableVariant(
                    var_id=principal.coordinates.string_format,
                )
            )

        return classifications
