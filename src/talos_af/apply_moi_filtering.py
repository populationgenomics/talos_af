"""
Ok, this is the stripped down Talos'y section

 - parse the VCF of labelled variants
 - group the variants by gene (for genes where recessive inheritance is relevant)
 - process each gene in turn according to its own MOI in the parsed specification
"""

import json
from argparse import ArgumentParser
from collections import defaultdict

from mendelbrot.pedigree_parser import PedigreeParser

from talos_af import check_moi, models
from talos_af import utils as utils_af


def set_up_filters(parsed_spec: dict, pedigree: PedigreeParser) -> dict[str, check_moi.BaseMoi]:
    """Create one instance of each MOI interpretation class."""
    moi_dict: dict[str, check_moi.BaseMoi] = {}

    for gene_details in parsed_spec.values():
        if (moi := gene_details.get('moi')) not in moi_dict:
            match moi:
                case 'AD':
                    moi_dict[moi] = check_moi.AD(pedigree)
                case 'AR':
                    moi_dict[moi] = check_moi.AR(pedigree)
                case 'XL':
                    moi_dict[moi] = check_moi.XL(pedigree)
    return moi_dict


def main(
    vcf_path: str,
    acmg_spec_path: str,
    pedigree_path: str,
    output_path: str,
    prior: str | None = None,
):
    """All the things."""

    with open(acmg_spec_path) as f:
        acmg_spec = json.load(f)

    id_lookup = {value['gene']: key for key, value in acmg_spec.items()}

    pedigree = PedigreeParser(pedigree_path)

    moi_filter_dict = set_up_filters(acmg_spec, pedigree=pedigree)

    # gather all variants indexed by gene
    gene_dict = utils_af.gather_gene_dict_from_vcf(vcf_path, id_lookup)

    selected_variants: dict[str, models.VariantAf] = {}
    sample_results: dict[str, list[models.ReportableVariant]] = defaultdict(list)

    for gene_id, variants in gene_dict.items():
        comp_het_dict = utils_af.find_comp_hets(variants, pedigree)
        moi_to_use = acmg_spec[gene_id]['moi']
        for variant in variants:
            if results := moi_filter_dict[moi_to_use].run(variant, comp_het_dict):
                # apply specific rules if required
                if acmg_spec[gene_id]['reportable'] == 'specific':
                    satisfies_conditions = utils_af.apply_gene_specific_rules(
                        rule=acmg_spec[gene_id]['specific_type'],
                        variant=variant,
                    )
                    if not satisfies_conditions:
                        continue

                selected_variants[variant.coordinates.string_format] = variant
                for sample, instances in results.items():
                    sample_results[sample].extend(instances)

    results_from_this_run = models.ResultsAf(variants=selected_variants, instances=sample_results)

    # if there are prior results, use those to reset first-seen dates
    utils_af.update_dates_from_prior_data(results_from_this_run, prior_path=prior)

    # write the output to long term storage using Pydantic - write if data is valid against the schema, write to file
    with open(output_path, 'w') as out_file:
        out_file.write(models.ResultsAf.model_validate(results_from_this_run).model_dump_json(indent=4))


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--vcf', help='Labelled and annotated VCF file', required=True)
    parser.add_argument('--acmg_spec', help="Specification of each gene's interpretation rules", required=True)
    parser.add_argument('--pedigree', help='Pedigree for the callset (required for sex)', required=True)
    parser.add_argument('--output', help='Path to write the output results to', required=True)
    parser.add_argument('--prior_results', help='Optional, Path to a previous set of results')
    args = parser.parse_args()
    main(
        vcf_path=args.vcf,
        acmg_spec_path=args.acmg_spec,
        pedigree_path=args.pedigree,
        output_path=args.output,
        prior=args.prior_results,
    )
