"""
takes an alphamissense tsv as input, filters regions, and outputs a new... VCF file.

input columns:
#CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class

Using https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz?download=1
- 613MB, containing all the pre-computed data we're interested in
"""

import gzip
from argparse import ArgumentParser
from importlib import resources

from talos_af.utils import REGION_DICT, process_bed, region_of_interest


def parse_and_filter_tsv(input_file: str, regions: REGION_DICT, output: str):
    with (
        gzip.open(input_file, 'rt') as handle,
        gzip.open(output, 'wt') as out,
        resources.open_text('talos_af', 'am_header.txt') as head_in,
    ):
        for line in head_in:
            out.write(line)

        for line in handle:
            if line.startswith('#'):
                continue

            llist = line.rstrip().split()

            chrom = llist[0]
            position = llist[1]

            if region_of_interest(regions=regions, chrom=chrom, pos=int(position)):
                out.write(
                    f'{chrom}\t{position}\t.\t{llist[2]}\t{llist[3]}\t60\tPASS\tam_class={llist[9]};am_score={llist[8]};am_transcript={llist[6].split(".")[0]}\n'
                )


def main(input_am: str, regions: str, output: str):
    bed_lookup: REGION_DICT = process_bed(bed_file=regions)
    parse_and_filter_tsv(input_file=input_am, regions=bed_lookup, output=output)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='input gzipped alphamissense tsv')
    parser.add_argument('--regions', help='input bed file defining regions of interest')
    parser.add_argument('--output', help='output VCF')
    args = parser.parse_args()
    main(input_am=args.input, regions=args.regions, output=args.output)
