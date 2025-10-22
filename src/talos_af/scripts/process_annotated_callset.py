#!/usr/bin/env python3

"""
Takes the full un-annotated joint-callset VCF,
A HailTable of the formatted annotations,
Integrates the two, writing a MatrixTable representation of the fully annotated VCF
"""

import json
from argparse import ArgumentParser

import hail as hl
from loguru import logger

MISSING_FLOAT_LO = hl.float64(0.0)
MISSING_INT = hl.int32(0)
MISSING_STRING = hl.str('missing')
ONE_INT = hl.int32(1)
BENIGN = hl.str('benign')
PATHOGENIC = hl.str('Pathogenic/Likely Pathogenic')
CRITICAL_CSQ_DEFAULT = [
    'frameshift',
    'splice_acceptor',
    'splice_donor',
    'start_lost',
    'stop_gained',
    'stop_lost',
    'transcript_ablation',
]
CSQ_STRING = [
    'consequence',
    'gene_id',
    'gene',
    'transcript',
    'mane_id',
    'mane',
    'biotype',
    'dna_change',
    'amino_acid_change',
    'codon',
    'ensp',
    'am_class',
    'am_pathogenicity',
    'revel_score',
]


def annotate_clinvar(mt: hl.MatrixTable, clinvar: str) -> hl.MatrixTable:
    """
    Don't allow these annotations to be missing
    - Talos has been co-developed with ClinvArbitration, a ClinVar re-summary effort
    - We no longer permit this to be missing (this has slipped in the past, causing odd results)

    See: https://github.com/populationgenomics/ClinvArbitration

    We replace any existing ClinVar annotations with our own version, then identify Pathogenic/Benign variants
    """

    logger.info(f'loading private clinvar annotations from {clinvar}')
    ht = hl.read_table(clinvar)
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            clinvar_significance=hl.or_else(ht[mt.row_key].clinical_significance, MISSING_STRING),
            clinvar_stars=hl.or_else(ht[mt.row_key].gold_stars, MISSING_INT),
            clinvar_allele=hl.or_else(ht[mt.row_key].allele_id, MISSING_INT),
        ),
    )

    # remove all confidently benign
    mt = mt.filter_rows(
        (mt.info.clinvar_significance.lower().contains(BENIGN)) & (mt.info.clinvar_stars > 0),
        keep=False,
    )

    # annotate as either strong or regular, return the result
    return mt.annotate_rows(
        info=mt.info.annotate(
            clinvar_plp=hl.if_else(
                (mt.info.clinvar_significance == PATHOGENIC) & (mt.info.clinvar_stars > 0),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def split_rows_by_gene_and_filter_to_green(mt: hl.MatrixTable, acmg_genes: hl.SetExpression) -> hl.MatrixTable:
    """
    splits each GeneId onto a new row, then filters any rows not annotating a relevant AF gene

    - first explode the matrix, separate gene per row
    - throw away all rows without a relevant gene
    - on all remaining rows, filter transcript consequences to match _this_ gene
    """

    # split each gene onto a separate row, transforms 'gene_ids' field from set to string
    mt = mt.explode_rows(mt.gene_ids)

    # filter rows without a green gene (removes empty gene_ids)
    mt = mt.filter_rows(acmg_genes.contains(mt.gene_ids))

    # limit the per-row transcript CSQ to those relevant to the single
    # gene now present on each row
    return mt.annotate_rows(
        transcript_consequences=mt.transcript_consequences.filter(lambda x: (mt.gene_ids == x.gene_id)),
    )


def annotate_category_high_impact(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    applies the boolean highimpact flag
    - Critical protein consequence on at least one transcript
    """

    critical_consequences = hl.set(CRITICAL_CSQ_DEFAULT)

    # First check if we have any HIGH consequences
    return mt.annotate_rows(
        info=mt.info.annotate(
            highimpact=hl.if_else(
                (
                    hl.len(
                        mt.transcript_consequences.filter(
                            lambda x: (
                                hl.len(critical_consequences.intersection(hl.set(x.consequence.split('&')))) > 0
                            ),
                        ),
                    )
                    > 0
                ),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def filter_to_population_rare(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    run the rare filter, using Gnomad Exomes and Genomes
    allow clinvar pathogenic to slip through this filter
    """
    rare_af_threshold = hl.literal(0.01)
    return mt.filter_rows(
        (mt.gnomad.gnomad_AF < rare_af_threshold) | (mt.info.clinvar_plp == ONE_INT),
    )


def csq_struct_to_string(tx_expr: hl.expr.StructExpression) -> hl.expr.ArrayExpression:
    """
    Taken shamelessly from the gnomad library source code
    Given a CSQ Struct, returns an array of CSQ strings
    (1 per csq in the struct).
    Fields & order correspond to those in `csq_fields`, corresponding to the
    VCF header that is required to interpret the VCF CSQ INFO field.
    Order is flexible & all fields in the default value are supported.
    These fields are formatted in the same way that their VEP CSQ counterparts are.

    Args:
        tx_expr (hl.Struct):
    Returns:
        generates an array of Strings for each CSQ
    """

    def get_csq_from_struct(element: hl.expr.StructExpression) -> hl.expr.StringExpression:
        # Most fields are 1-1, just lowercase
        fields = dict(element)
        return hl.delimit([hl.or_else(hl.str(fields.get(f, '')), '') for f in CSQ_STRING], '|')

    csq = hl.empty_array(hl.tstr)
    csq = csq.extend(
        hl.or_else(tx_expr.map(lambda x: get_csq_from_struct(x)), hl.empty_array(hl.tstr)),
    )

    # previous consequence filters may make this caution unnecessary
    return hl.or_missing(hl.len(csq) > 0, csq)


def write_matrix_to_vcf(mt: hl.MatrixTable, vcf_out: str):
    """
    write the remaining MatrixTable content to file as a VCF
    generate a custom header containing the CSQ contents which
    were retained during this run

    Args:
        mt (): the whole MatrixTable
        vcf_out (str): where to write the VCF
    """

    # this temp file needs to be in GCP, not local
    # otherwise the batch that generates the file won't be able to read
    header_path = 'additional_header.txt'

    # generate a CSQ string specific to the config file for decoding later
    csq_contents = '|'.join(CSQ_STRING)

    # write this custom header locally
    with open(header_path, 'w') as handle:
        handle.write(f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {csq_contents}">')
    logger.info(f'Writing categorised variants out to {vcf_out}')
    hl.export_vcf(mt, vcf_out, append_to_header=header_path, tabix=True)


def cli_main():
    """
    take an input and output VCF path -
    - input path contains all variants
    - this is combined with a HT of annotations to make an annotated callset
    - the annotations are applied to generate a full MatrixTable of annotations and variants
    - this is then filtered down to plausibly pathogenic variants

    Does this need to be in Hail? Probably not, but it lets me re-use as much logic as possible
    """

    parser = ArgumentParser(description='Takes a HT of annotations, and a callset VCF, and combines into a MT')
    parser.add_argument('--input', help='Path to the MatrixTable', required=True)
    parser.add_argument('--annotations', help='Path to the annotated sites-only VCF', required=True)
    parser.add_argument('--clinvar', help='Path to the ClinVar HT', required=True)
    parser.add_argument('--acmg_spec', help='Path to the ACMG spec JSON', required=True)
    parser.add_argument('--output', help='output Table path, must have a ".ht" extension', required=True)
    args = parser.parse_args()

    logger.info(
        r"""Welcome To
    _________ _______  _        _______  _______         _______  _______
    \__   __/(  ___  )( \      (  ___  )(  ____ \       (  ___  )(  ____ \
       ) (   | (   ) || (      | (   ) || (    \/       | (   ) || (    \/
       | |   | (___) || |      | |   | || (_____  _____ | (___) || (__
       | |   |  ___  || |      | |   | |(_____  )(_____)|  ___  ||  __)
       | |   | (   ) || |      | |   | |      ) |       | (   ) || (
       | |   | )   ( || (____/\| (___) |/\____) |       | )   ( || )
       )_(   |/     \|(_______/(_______)\_______)       |/     \||/
   """,
    )

    main(
        input_path=args.input,
        annotations=args.annotations,
        clinvar_path=args.clinvar,
        acmg_spec_path=args.acmg_spec,
        output_path=args.output,
    )


def main(
    input_path: str,
    annotations: str,
    clinvar_path: str,
    acmg_spec_path: str,
    output_path: str,
):
    """
    Takes a Hail-Table of annotations, a joint-called VCF, reads the VCF as a MatrixTable and hops the annotations over

    Args:
        input_path (str): path to the full callset VCF, or MT of the same
        annotations (str): path to a Hail Table containing annotations
        clinvar_path (str): path to a Hail Table containing ClinVar annotations
        acmg_spec_path (str): path to a JSON file containing ACMG annotations
        output_path (str): path to write the resulting MatrixTable to
    """

    logger.info('Using local backend for Hail')
    hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)

    # read the VCF into a MatrixTable, or read the existing MatrixTable
    if input_path.endswith('.mt'):
        logger.info(f'Reading MatrixTable from {input_path!r}')
        mt = hl.read_matrix_table(input_path)
    elif input_path.endswith(('.vcf', '.vcf.gz', '.vcf.bgz')):
        logger.info(f'Reading VCF from {input_path!r}')
        # import the VCF, and set the array_elements_required to False
        # this is because the VCF may not have all fields present in all rows
        # which is fine for our purposes
        mt = hl.import_vcf(
            input_path,
            array_elements_required=False,
            force_bgz=True,
        )
    else:
        raise ValueError(f'Input path must be a VCF or MatrixTable, got {input_path!r}')

    with open(acmg_spec_path) as f:
        acmg_spec = json.load(f)
    acmg_genes = hl.literal(set(acmg_spec.keys()))

    # read the annotations into a Table
    ht = hl.read_table(annotations)

    # syntax sweeter for later on
    matched_annotations = ht[mt.row_key]

    # a couple of lines commented of to make this as easy as possible to adopt
    # talos doesn't make use of these annotations yet
    mt = mt.annotate_rows(
        gnomad=hl.struct(
            gnomad_AC=matched_annotations.info.gnomad_AC_joint,
            gnomad_AF=matched_annotations.info.gnomad_AF_joint,
            gnomad_AC_XY=matched_annotations.info.gnomad_AC_joint_XY,
            gnomad_HomAlt=matched_annotations.info.gnomad_HomAlt_joint,
        ),
        transcript_consequences=matched_annotations.transcript_consequences,
        gene_ids=matched_annotations.gene_ids,
    )

    mt = annotate_clinvar(mt, clinvar_path)

    # checkpoint the MatrixTable locally to make everything downstream faster
    mt = mt.checkpoint('checkpoint.mt', overwrite=True, _read_if_exists=True)

    # now do some seletion/minimisation
    mt = split_rows_by_gene_and_filter_to_green(mt, acmg_genes)

    # find high impact genes
    mt = annotate_category_high_impact(mt)

    # and now filter down to relevant variants, either being ClinVar pathogenic or high impact
    mt = mt.filter_rows((mt.info.clinvar_plp == 1) | (mt.info.highimpact == 1))

    # obtain the massive CSQ string using method stolen from the Broad's Gnomad library
    # also take the single gene_id (from the exploded attribute)
    # and retrieves the gnomAD annotations
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            **mt.gnomad,
            csq=csq_struct_to_string(mt.transcript_consequences),
            gene_id=mt.gene_ids,
        ),
    )

    write_matrix_to_vcf(mt=mt, vcf_out=output_path)


if __name__ == '__main__':
    cli_main()
