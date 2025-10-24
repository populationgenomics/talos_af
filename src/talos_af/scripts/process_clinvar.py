"""
takes a clinvar csv as input, filters regions, and outputs a new... BED file. TBC

input columns:
contig  position        reference       alternate       clinical_significance   gold_stars      allele_id
"""

import gzip
from argparse import ArgumentParser
from csv import DictReader

from talos_af.utils import REGION_DICT, process_bed, region_of_interest


def parse_and_filter_tsv(input_file: str, regions: REGION_DICT, header: str):
    """Read the compressed TSV input file, filter it by the acceptable regions, and write as a new minimised BED."""
    with open(input_file) as handle, open(header, 'r') as head_in:
        for line in head_in:
            print(line.rstrip())

        for line in DictReader(handle, delimiter='\t'):
            chrom = line['contig']
            pos = line['position']
            if (
                region_of_interest(regions=regions, chrom=chrom, pos=int(pos))
                and (clinsig := line['clinical_significance']) != 'VUS'
            ):
                ref = line['reference']
                alt = line['alternate']
                sig = f'clinical_significance={clinsig}'
                stars = f'gold_stars={line["gold_stars"]}'
                allele_id = f'allele_id={line["allele_id"]}'
                print(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t{sig};{stars};{allele_id}')


def main(input: str, regions: str, header: str):
    bed_lookup = process_bed(bed_file=regions)
    parse_and_filter_tsv(input_file=input, regions=bed_lookup, header=header)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='input clinvar tsv')
    parser.add_argument('--regions', help='input bed file containing regions of interest')
    parser.add_argument('--header', help='header file for the output VCF')
    args = parser.parse_args()
    main(input=args.input, regions=args.regions, header=args.header)
