from collections import defaultdict
from itertools import combinations_with_replacement

import cyvcf2
from loguru import logger
from mendelbrot.pedigree_parser import PedigreeParser

from talos_af import models
from talos_af.scripts.process_annotated_callset import CSQ_STRING

REGION_DICT = dict[str, list[tuple[int, int]]]

HOMREF: int = 0
HETALT: int = 1
UNKNOWN: int = 2
HOMALT: int = 3

# in cyVCF2, these ints represent HOMREF, and UNKNOWN
BAD_GENOTYPES: set[int] = {HOMREF, UNKNOWN}
NON_HOM_CHROM = ['X', 'Y', 'MT', 'M']
X_CHROMOSOME = {'X'}

# used by cyvcf2 for... not phased
PHASE_SET_DEFAULT = -2147483648

CATEGORIES = {'highimpact', 'clinvar_plp'}
GeneDict = dict[str, list[models.VariantAf]]
CompHetDict = dict[str, dict[str, list[models.VariantAf]]]


def process_bed(bed_file: str) -> REGION_DICT:
    """
    Read a BED file and process it into a series of intervals

    Args:
        bed_file: str

    Returns:
        a hierarchically organised dictionary of regions
    """

    bed_regions: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(bed_file, encoding='utf-8') as handle:
        for line in handle:
            llist = line.strip().split('\t')
            chrom = llist[0]
            start = int(llist[1])
            end = int(llist[2])
            bed_regions[chrom].append((start, end))

    return bed_regions


def region_of_interest(regions: REGION_DICT, chrom: str, pos: int) -> bool:
    """Check if this is a region of interest."""
    if chrom not in regions:
        return False

    for start, end in regions[chrom]:
        if start <= pos <= end:
            return True
    return False


def extract_csq(csq_contents: str) -> list[dict]:
    """Handle extraction of the CSQ entries."""

    # allow for no CSQ data, i.e. splice variant
    if not csq_contents:
        return []

    # iterate over all consequences, and make each into a dict
    txc_dicts = [dict(zip(CSQ_STRING, each_csq.split('|'), strict=True)) for each_csq in csq_contents.split(',')]

    # update this String to be either a float, or missing
    for each_dict in txc_dicts:
        am_path = each_dict.get('am_pathogenicity')
        each_dict['am_pathogenicity'] = float(am_path) if am_path else ''

    return txc_dicts


def get_non_ref_samples(variant: 'cyvcf2.Variant', samples: list[str]) -> tuple[set[str], set[str]]:
    """
    For this variant, find all samples with a call. cyvcf2 uses 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT.
    Returns 2 sets of strings; het sample ID, hom sample IDs
    """
    het_samples = set()
    hom_samples = set()

    # this iteration is based on the cyvcf2 representations
    for sam, genotype_int in zip(samples, variant.gt_types, strict=True):
        if genotype_int in BAD_GENOTYPES:
            continue
        if genotype_int == HETALT:
            het_samples.add(sam)
        if genotype_int == HOMALT:
            hom_samples.add(sam)

    return het_samples, hom_samples


def get_phase_data(samples: list[str], var: 'cyvcf2.Variant') -> dict[str, dict[int, str]]:
    """
    read phase data from this variant

    we tolerate a variety of failures, the worst case scenario is that we have no phase data
    PS - PhaseSet
    PID/PGT - PhaseID/PhaseGenotype

    No other phase types are currently supported

    Args:
        samples (list[str]): all samples in the VCF
        var (cyvcf2.Variant):
    """

    global PHASE_BROKEN

    phased_dict: dict[str, dict[int, str]] = defaultdict(dict)

    # check we have relevant attributes in the variant format fields
    if not any(x in var.FORMAT for x in ['PS', 'PGT', 'PID']):
        return dict(phased_dict)

    # first set the numpy.ndarray to be a list of ints
    # then zip against ordered sample IDs
    # this might need to store the exact genotype too
    # i.e. 0|1 and 1|0 can be in the same phase-set
    # but are un-phased variants
    try:
        if 'PS' in var.FORMAT:
            for sample, phase, genotype in zip(samples, map(int, var.format('PS')), var.genotypes, strict=True):
                # cyvcf2.Variant holds two ints, and a bool for biallelic calls
                # but only one int and a bool for hemi
                if len(genotype) == 3:
                    allele_1, allele_2, phased = genotype
                    gt = f'{allele_1}|{allele_2}'
                elif len(genotype) == 2:
                    allele_1, phased = genotype
                    gt = f'{allele_1}'
                else:
                    raise ValueError(f'Unexpected genotype length: {len(genotype)}: {genotype}')

                if not phased:
                    continue

                # phase set is a number
                if phase != PHASE_SET_DEFAULT:
                    phased_dict[sample][phase] = gt

        elif all(attr_name in var.FORMAT for attr_name in ['PGT', 'PID']):
            logger.info('Failed to find PS phase attributes')
            try:
                # retry using PGT & PID
                for sample, phase_gt, phase_id in zip(
                    samples,
                    var.format('PGT'),
                    var.format('PID'),
                    strict=True,
                ):
                    if phase_gt != '.' and phase_id != '.':
                        phased_dict[sample][phase_id] = phase_gt

            except KeyError as ke2:
                logger.info('Failed to determine phase information using PID and PGT')
                raise ke2
        elif not PHASE_BROKEN:
            logger.info('Found no PS phase attributes (known formats are PS, PGT/PID)')
            PHASE_BROKEN = True

    except (KeyError, ValueError):
        if not PHASE_BROKEN:
            logger.info('Failed to correctly parse known phase attributes using existing methods')
            logger.info(
                'Please post an issue on the Talos GitHub Repo with the VCF FORMAT lines and descriptions',
            )
        PHASE_BROKEN = True

    return dict(phased_dict)


def create_small_variant(
    var: 'cyvcf2.Variant',
    samples: list[str],
) -> models.VariantAf | None:
    """Takes a small variant and creates a Model from it."""

    coordinates = models.Coordinates(chrom=var.CHROM.replace('chr', ''), pos=var.POS, ref=var.REF, alt=var.ALT[0])
    info: dict[str, str | int | float] = {x.lower(): y for x, y in var.INFO} | {'var_link': coordinates.string_format}

    het_samples, hom_samples = get_non_ref_samples(variant=var, samples=samples)

    # overwrite with true booleans
    for cat in CATEGORIES:
        info[cat] = info.get(cat, 0) == 1

    phased = get_phase_data(samples, var)

    # only keep these where the sample has a variant - the majority of samples have empty data, and we don't use it
    # if we require depths/ratios/etc. for WT samples, revisit this
    # hopefully a solution to the memory explosion in large cohorts
    transcript_consequences = extract_csq(csq_contents=info.pop('csq', ''))

    return models.VariantAf(
        coordinates=coordinates,
        info=info,
        het_samples=het_samples,
        hom_samples=hom_samples,
        phased=phased,
        transcript_consequences=transcript_consequences,
    )


def gather_gene_dict_from_vcf(vcf_path: str) -> dict[str, list[models.VariantAf]]:
    """
    takes a cyvcf2.VCFReader instance, and a specified chromosome
    iterates over all variants in the region, and builds a lookup

    optionally takes a second VCF and incorporates into same dict

    Args:
        vcf_path (str): the VCF to read

    Returns:
        A lookup in the form
        {
            gene1: [var1, var2],
            gene2: [var3],
            ...
        }
    """

    # a dict to allow lookup of variants on this whole chromosome
    variant_count = 0
    contig_dict = defaultdict(list)

    variant_source = cyvcf2.VCF(vcf_path)

    # iterate over all variants on this contig and store by unique key
    # if contig has no variants, prints an error and returns []
    for variant_row in variant_source:
        variant = create_small_variant(var=variant_row, samples=variant_source.samples)
        # update the variant count
        variant_count += 1
        # update the gene index dictionary
        contig_dict[variant.info['gene_id']].append(variant)

    logger.info(f'VCF {vcf_path} contained {variant_count} variants, in {len(contig_dict)} genes')

    return dict(contig_dict)


def find_comp_hets(var_list: list[models.VariantAf], pedigree: PedigreeParser) -> CompHetDict:
    """
    Find compound het pairs, variants provided in the format [var1, var2, ...]

    generates pair content in the form
    {
        sample: {
            var_as_string: [partner_variant, ...],
        }
    }
    """

    # create an empty dictionary
    comp_het_results: CompHetDict = defaultdict(dict)

    # use combinations_with_replacement to find all gene pairs
    for var_1, var_2 in combinations_with_replacement(var_list, 2):
        assert var_1.coordinates.chrom == var_2.coordinates.chrom

        if (var_1.coordinates == var_2.coordinates) or var_1.coordinates.chrom in NON_HOM_CHROM:
            continue

        # iterate over any samples with a het overlap
        for sample in var_1.het_samples.intersection(var_2.het_samples):
            phased = False

            # don't assess male compound hets on sex chromosomes
            if pedigree.participants[sample].sex == 1 and var_1.coordinates.chrom in X_CHROMOSOME:
                continue

            # check for both variants being in the same phase set
            if sample in var_1.phased and sample in var_2.phased:
                # check for presence of the same phase set
                for phase_set in [ps for ps in var_1.phased[sample] if ps in var_2.phased[sample]]:
                    if var_1.phased[sample][phase_set] == var_2.phased[sample][phase_set]:
                        phased = True
            if not phased:
                comp_het_results[sample].setdefault(var_1.coordinates.string_format, []).append(var_2)
                comp_het_results[sample].setdefault(var_2.coordinates.string_format, []).append(var_1)

    return comp_het_results
