"""
takes a revel csv as input, filters regions, and outputs a new... BED file. TBC

input columns:
chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid

trying this with https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip?download=1
- this is the 630MB file, compressing 6.5GB of raw pre-computed data
"""

import gzip
import zipfile
from argparse import ArgumentParser

from talos_af.utils import REGION_DICT, process_bed, region_of_interest


def parse_and_filter_tsv(input_file: str, regions: REGION_DICT, header: str, output: str):
    """Read the compressed TSV input file, filter it by the acceptable regions, and write as a new minimised BED."""

    with zipfile.ZipFile(input_file) as ziphandle, gzip.open(output, 'wt') as out, open(header, 'r') as head_in:
        for line in head_in:
            out.write(line)

        filename = ziphandle.namelist()[0]
        with ziphandle.open(filename, 'r') as handle:
            for line in map(str, handle.__iter__()):
                # there simply must be a neater way of doing this, it's horrible
                trimline = line.lstrip("b'").rstrip("\\n'")
                llist = trimline.split(',')

                chrom = 'chr' + llist[0]
                position = llist[2]

                # there are some of these
                if position == '.' or chrom == 'chrchr':
                    continue

                if region_of_interest(regions=regions, chrom=chrom, pos=int(position)):
                    # in the TSV these are semi-colon separated which wouldn't be permitted in the VCF
                    transcripts = llist[8].replace(';', ',')
                    out.write(
                        f'{chrom}\t{position}\t.\t{llist[3]}\t{llist[4]}\t60\tPASS\tREVEL={llist[7]}~{transcripts}\n'
                    )


def main(input_revel: str, regions: str, header: str, output: str):
    bed_lookup = process_bed(bed_file=regions)
    parse_and_filter_tsv(input_file=input_revel, regions=bed_lookup, header=header, output=output)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='input compressed revel csv')
    parser.add_argument('--regions', help='input bed file containing regions of interest')
    parser.add_argument('--header', help='header file for the output VCF')
    parser.add_argument('--output', help='output path for the resulting VCF file')
    args = parser.parse_args()
    main(input_revel=args.input, regions=args.regions, header=args.header, output=args.output)
