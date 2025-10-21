#!/usr/bin/env python3

"""
This is an adapter process to take a sites-only VCF annotated with gnomAD frequencies and BCFtools csq consequences, and
re-arrange it into a HailTable for use with the Talos pipeline.

This process combines the AF/CSQs already applied with REVEL and AlphaMissense annotations
"""

import json
from argparse import ArgumentParser

import hail as hl
from loguru import logger

MISSING_FLOAT = hl.float64(0)
MISSING_INT = hl.int32(0)
MISSING_STRING = hl.str('')


def extract_and_split_csq_string(vcf_path: str) -> list[str]:
    """
    Extract the BCSQ header from the VCF and split it into a list of strings

    Args:
        vcf_path (str): path to the local VCF

    Returns:
        list of strings
    """

    # get the headers from the VCF
    all_headers = hl.get_vcf_metadata(vcf_path)

    # get the '|'-delimited String of all header names
    csq_whole_string = all_headers['info']['BCSQ']['Description'].split('Format: ')[-1]

    # split it all on pipes, return the list
    return csq_whole_string.lower().split('|')


def csq_strings_into_hail_structs(csq_strings: list[str], ht: hl.Table | hl.MatrixTable) -> hl.Table | hl.MatrixTable:
    """
    Take the list of BCSQ strings, split the CSQ annotation and re-organise as a hl struct

    Args:
        csq_strings (list[str]): a list of strings, each representing a CSQ entry
        ht (hl.Table): the Table to annotate

    Returns:
        a Table with the BCSQ annotations re-arranged
    """

    # get the BCSQ contents as a list of lists of strings, per variant
    split_csqs = ht.info.BCSQ.map(lambda csq_entry: csq_entry.split('\|'))  # noqa: W605

    # this looks pretty hideous, bear with me
    # if BCFtools csq doesn't have a consequence annotation, it will truncate the pipe-delimited string
    # this is fine sometimes, but not when we're building a schema here
    # when we find truncated BCSQ strings, we need to add dummy values to the end of the array
    split_csqs = split_csqs.map(
        lambda x: hl.if_else(
            # if there were only 4 values, add 3 missing Strings
            hl.len(x) == 4,
            x.extend([MISSING_STRING, MISSING_STRING, MISSING_STRING]),
            hl.if_else(
                # 5 values... add 2 missing Strings
                hl.len(x) == 5,
                x.extend([MISSING_STRING, MISSING_STRING]),
                hl.if_else(
                    hl.len(x) == 6,
                    x.extend([MISSING_STRING]),
                    x,
                ),
            ),
        ),
    )

    # transform the CSQ string arrays into structs using the header names
    # Consequence | gene | transcript | biotype | strand | amino_acid_change | dna_change
    ht = ht.annotate(
        transcript_consequences=split_csqs.map(
            lambda x: hl.struct(
                **{csq_strings[n]: x[n] for n in range(len(csq_strings)) if csq_strings[n] != 'strand'},
            ),
        ),
    )

    return ht.annotate(
        # amino_acid_change can be absent, or in the form of "123P" or "123P-124F"
        # we use this number when matching to the codons of missense variants, to find codon of the reference pos.
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                codon=hl.if_else(
                    x.amino_acid_change == MISSING_STRING,
                    hl.missing(hl.tint32),
                    hl.if_else(
                        x.amino_acid_change.matches('^([0-9]+).*$'),
                        hl.int32(x.amino_acid_change.replace('^([0-9]+).+', '$1')),
                        hl.missing(hl.tint32),
                    ),
                ),
            ),
            ht.transcript_consequences,
        ),
    )


def annotate_gene_ids(ht: hl.Table, acmg_spec_path: str) -> hl.Table:
    """
    Using a JSON file as a lookup, annotate transcript consequences with Ensembl gene ID.
    """

    with open(acmg_spec_path) as f:
        acmg_spec = json.load(f)

    id_hl_dict = hl.literal({value['gene_id']: key for key, value in acmg_spec.items()})

    # take the ENSG value from the dict for the contig (correctly matches PAR region genes)
    # default to the gene symbol (which can be the ENSG, depending on transcript consequence)
    return ht.annotate(
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                gene_id=id_hl_dict.get(x.gene, x.gene),
            ),
            ht.transcript_consequences,
        ),
    )


def insert_ext_annotations(ht: hl.Table, am_table_path: str, revel_table_path: str) -> hl.Table:
    """
    Load up Hail Tables of External annotations, add this data by matching transcript IDs.
    """

    logger.info(f'Reading AM annotations from {am_table_path} and applying to MT')
    logger.info(f'Reading REVEL annotations from {revel_table_path} and applying to MT')

    # read in the hail table containing alpha missense annotations
    am_ht = hl.read_table(am_table_path)
    revel_ht = hl.read_table(revel_table_path)

    # AM consequence matching needs conditional application based on the specific transcript match
    return ht.annotate(
        transcript_consequences=hl.map(
            lambda x: x.annotate(
                am_class=hl.if_else(
                    x.transcript == am_ht[ht.key].transcript,
                    am_ht[ht.key].am_class,
                    MISSING_STRING,
                ),
                am_score=hl.if_else(
                    x.transcript == am_ht[ht.key].transcript,
                    am_ht[ht.key].am_score,
                    MISSING_FLOAT,
                ),
                revel_score=hl.if_else(
                    x.transcript == revel_ht[ht.key].transcript,
                    revel_ht[ht.key].score,
                    MISSING_STRING,
                ),
            ),
            ht.transcript_consequences,
        ),
    )


def cli_main():
    """
    take an input VCF and an output MT path
    also supply the alpha_missense table created by parse_amissense_into_ht.py
    """

    parser = ArgumentParser(description='Takes a BCSQ annotated VCF and makes it a HT')
    parser.add_argument('--input', help='Path to the annotated sites-only VCF', required=True)
    parser.add_argument('--acmg_spec', help='BED file containing gene mapping')
    parser.add_argument('--am', help='Hail Table containing AlphaMissense annotations', required=True)
    parser.add_argument('--revel', help='Hail Table containing REVEL annotations', required=True)
    parser.add_argument('--output', help='output Table path, must have a ".ht" extension', required=True)
    args = parser.parse_args()

    main(
        vcf_path=args.input,
        output_path=args.output,
        acmg_spec=args.acmg_spec,
        alpha_m=args.am,
        revel=args.revel,
    )


def main(
    vcf_path: str,
    output_path: str,
    acmg_spec: str,
    alpha_m: str,
    revel: str,
):
    """
    Takes a BCFtools-annotated VCF, reorganises into a Talos-compatible MatrixTable
    Will annotate at runtime with AlphaMissense annotations

    Args:
        vcf_path (str): path to the annotated sites-only VCF
        output_path (str): path to write the resulting Hail Table to, must
        acmg_spec (str): path to an acmg_spec JSON
        alpha_m (str): path to the AlphaMissense Hail Table, required
        revel (str): path to a REVEL Hail Table for enhanced annotation
    """

    logger.info('Using Hail Local backend')
    hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)

    # pull and split the CSQ header line
    csq_fields = extract_and_split_csq_string(vcf_path=vcf_path)

    # read the VCF into a MatrixTable
    mt = hl.import_vcf(vcf_path, array_elements_required=False, force_bgz=True)

    # checkpoint the rows as a Table locally to make everything downstream faster
    ht = mt.rows().checkpoint('checkpoint.ht', overwrite=True, _read_if_exists=True)

    logger.info('VCF imported and checkpointed as a Hail Table')

    # re-shuffle the BCSQ elements
    ht = csq_strings_into_hail_structs(csq_fields, ht)

    # add ENSG IDs
    ht = annotate_gene_ids(ht, acmg_spec_path=acmg_spec)

    # get a hold of the geneIds - use some aggregation
    ht = ht.annotate(gene_ids=hl.set(ht.transcript_consequences.map(lambda c: c.gene_id)))

    # checkpoint before combining with external tables
    ht = mt.rows().checkpoint('checkpoint_ext_tables.ht', overwrite=True, _read_if_exists=True)
    logger.info('Checkpointed prior to AM/Revel annotation')

    # add AlphaMissense scores
    ht = insert_ext_annotations(ht, am_table_path=alpha_m, revel_table_path=revel)

    # drop the BCSQ field
    ht = ht.annotate(info=ht.info.drop('BCSQ'))

    ht.describe()

    ht.write(output_path, overwrite=True)


if __name__ == '__main__':
    cli_main()
