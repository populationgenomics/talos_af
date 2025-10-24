"""
takes a clinvar csv as input, filters regions, and outputs a new... BED file. TBC

input columns:
contig  position        reference       alternate       clinical_significance   gold_stars      allele_id
"""

from argparse import ArgumentParser
from csv import DictReader

from talos_af.utils import REGION_DICT, process_bed, region_of_interest


def parse_and_filter_tsv(input_file: str, regions: REGION_DICT, header: str):
    """Read the compressed TSV input file, filter it by the acceptable regions, and write as a new minimised BED."""
    with open(input_file) as handle, open(header) as head_in:
        for line in head_in:
            print(line.rstrip())

        for dict_line in DictReader(handle, delimiter='\t'):
            if not isinstance(dict_line, dict):
                continue
            chrom = dict_line['contig']
            pos = dict_line['position']
            if (
                region_of_interest(regions=regions, chrom=chrom, pos=int(pos))
                and (clinsig := dict_line['clinical_significance']) != 'VUS'
            ):
                ref = dict_line['reference']
                alt = dict_line['alternate']
                sig = f'clinical_significance={clinsig}'
                stars = f'gold_stars={dict_line["gold_stars"]}'
                allele_id = f'allele_id={dict_line["allele_id"]}'
                print(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t{sig};{stars};{allele_id}')


def main(input_path: str, regions: str, header: str):
    bed_lookup = process_bed(bed_file=regions)
    parse_and_filter_tsv(input_file=input_path, regions=bed_lookup, header=header)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='input clinvar tsv')
    parser.add_argument('--regions', help='input bed file containing regions of interest')
    parser.add_argument('--header', help='header file for the output VCF')
    args = parser.parse_args()
    main(input_path=args.input, regions=args.regions, header=args.header)
