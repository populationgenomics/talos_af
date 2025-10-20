"""
takes a alphamissense tsv as input, filters regions, and outputs a new... file. TBC

input columns:
#CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class
"""

import gzip
from argparse import ArgumentParser
from collections import defaultdict

REGION_DICT = dict[str, list[tuple[int, int]]]


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
            bed_regions[llist[0]].append((int(llist[1]), int(llist[2])))

    return bed_regions


def region_of_interest(regions: REGION_DICT, chrom: str, pos: int) -> bool:
    """Check if this is a region of interest."""
    if chrom not in regions:
        return False

    for start, end in regions[chrom]:
        return start <= pos <= end
    return False


def main(input_am: str, input_bed: str, output: str):
    bed_lookup = process_bed(bed_file=input_bed)

    with gzip.open(input_am, 'rt') as handle, open(output, 'w') as out:
        for line in handle:
            if line.startswith('#'):
                continue

            # there simply must be a neater way of doing this, it's horrible
            llist = line.rstrip().split()

            chrom = llist[0]
            position = llist[1]

            if region_of_interest(regions=bed_lookup, chrom=chrom, pos=int(position)):
                out.write(
                    '\t'.join(
                        [
                            chrom,
                            position,
                            llist[2],  # ref allele
                            llist[3],  # alt allele
                            llist[6].split('.')[0],  # transcript, trim off the version portion
                            llist[8],  # score
                            llist[9],  # class
                        ]
                    ) + '\n'
                )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='input gzipped alphamissense tsv')
    parser.add_argument('--bed', help='input bed file')
    parser.add_argument('--output', help='output file')
    args = parser.parse_args()
    main(input_am=args.input, input_bed=args.bed, output=args.output)
