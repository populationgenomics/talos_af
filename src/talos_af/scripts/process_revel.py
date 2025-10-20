"""
takes a revel csv as input, filters regions, and outputs a new... file. TBC

input columns:
chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
"""

import zipfile
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
            chrom = llist[0][3:]
            bed_regions[chrom].append((int(llist[1]), int(llist[2])))

    return bed_regions


def region_of_interest(regions: REGION_DICT, chrom: str, pos: int) -> bool:
    """Check if this is a region of interest."""
    if chrom not in regions:
        return False

    for start, end in regions[chrom]:
        return start <= pos <= end
    return False


def main(input_revel: str, input_bed: str, output: str):
    bed_lookup = process_bed(bed_file=input_bed)

    with zipfile.ZipFile(input_revel) as ziphandle, open(output, 'w') as out:
        filename = ziphandle.namelist()[0]
        with ziphandle.open(filename, 'r') as handle:
            for line in map(str, handle.__iter__()):
                # there simply must be a neater way of doing this, it's horrible
                trimline = line.lstrip("b'").rstrip("\\n'")
                llist = trimline.split(',')

                chrom = llist[0]
                position = llist[2]

                # there are some of these
                if position == '.' or chrom == 'chr':
                    continue

                if region_of_interest(regions=bed_lookup, chrom=chrom, pos=int(position)):
                    out.write(
                        '\t'.join(
                            [
                                chrom,
                                position,
                                position,
                                llist[3],
                                llist[4],
                                llist[7],
                                llist[8],
                            ]
                        ) + '\n'
                    )


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='input compressed revel csv')
    parser.add_argument('--bed', help='input bed file')
    parser.add_argument('--output', help='output file')
    args = parser.parse_args()
    main(input_revel=args.input, input_bed=args.bed, output=args.output)
