"""
takes a alphamissense tsv as input, filters regions, and outputs a new... file. TBC

input columns:
#CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class

Using https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz?download=1
- 613MB, containing all the pre-computed data we're interested in
"""

import gzip
from argparse import ArgumentParser

import hail as hl

from talos_af.utils import REGION_DICT, process_bed, region_of_interest


def convert_to_ht(input_file: str, output_path: str):
    # start up a local hail sesh
    hl.context.init_spark(master='local[*]', default_reference='GRCh38', quiet=True)

    # read the TSV, setting some non-string types
    types = {'score': hl.tfloat32}

    ht = hl.import_table(input_file, types=types, delimiter='\t', skip_blank_lines=True)

    # now rearrange some bits
    ht = ht.transmute(
        locus=hl.parse_locus(ht.locus),
        alleles=[ht.ref, ht.alt],
    )

    # and set a key
    ht = ht.key_by('locus', 'alleles')
    ht.write(output_path)


def parse_and_filter_tsv(input_file: str, regions: REGION_DICT, output_path: str):
    with gzip.open(input_file, 'rt') as handle, open(output_path, 'w') as out:
        # write the file header
        out.write('locus\tref\talt\ttranscript\tscore\tclass\n')

        for line in handle:
            if line.startswith('#'):
                continue

            # there simply must be a neater way of doing this, it's horrible
            llist = line.rstrip().split()

            chrom = llist[0]
            position = llist[1]

            if region_of_interest(regions=regions, chrom=chrom, pos=int(position)):
                out.write(
                    '\t'.join(
                        [
                            f'{chrom}:{position}',
                            llist[2],  # ref allele
                            llist[3],  # alt allele
                            llist[6].split('.')[0],  # transcript, trim off the version portion
                            llist[8],  # score
                            llist[9],  # class
                        ]
                    )
                    + '\n'
                )


def main(input_am: str, input_bed: str, output_tsv: str, output_ht: str):
    bed_lookup: REGION_DICT = process_bed(bed_file=input_bed)
    parse_and_filter_tsv(input_file=input_am, regions=bed_lookup, output_path=output_tsv)
    convert_to_ht(input_file=output_tsv, output_path=output_ht)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='input gzipped alphamissense tsv')
    parser.add_argument('--bed', help='input bed file')
    parser.add_argument('--output_tsv', help='output TSV file')
    parser.add_argument('--output_ht', help='output HT path')
    args = parser.parse_args()
    main(input_am=args.input, input_bed=args.bed, output_tsv=args.output_tsv, output_ht=args.output_ht)
