"""
takes a revel csv as input, filters regions, and outputs a new... file. TBC

input columns:
chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid

trying this with https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip?download=1
- this is the 630MB file, compressing 6.5GB of raw pre-computed data
"""

import zipfile
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


def parse_and_filter_tsv(input_file: str, regions: REGION_DICT, output_tsv: str):
    """Read the compressed TSV input file, filter it by the acceptable regions, and write as a new minimised TSV."""

    with zipfile.ZipFile(input_file) as ziphandle, open(output_tsv, 'w') as out:
        # write the file header
        out.write('locus\tref\talt\tscore\ttranscript\n')

        filename = ziphandle.namelist()[0]
        with ziphandle.open(filename, 'r') as handle:
            for line in map(str, handle.__iter__()):
                # there simply must be a neater way of doing this, it's horrible
                trimline = line.lstrip("b'").rstrip("\\n'")
                llist = trimline.split(',')

                chrom = 'chr' + llist[0]
                position = llist[2]

                # there are some of these
                if position == '.' or chrom == 'chr':
                    continue

                if region_of_interest(regions=regions, chrom=chrom, pos=int(position)):
                    out.write(
                        '\t'.join(
                            [
                                f'chr{chrom}:{position}',
                                llist[3],
                                llist[4],
                                llist[7],
                                llist[8],
                            ]
                        )
                        + '\n'
                    )


def main(input_revel: str, input_bed: str, output_tsv: str, output_ht: str):
    bed_lookup = process_bed(bed_file=input_bed)
    parse_and_filter_tsv(input_file=input_revel, regions=bed_lookup, output_tsv=output_tsv)
    convert_to_ht(input_file=output_tsv, output_path=output_ht)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='input compressed revel csv')
    parser.add_argument('--bed', help='input bed file')
    parser.add_argument('--output_tsv', help='output TSV file')
    parser.add_argument('--output_ht', help='output HT path')
    args = parser.parse_args()
    main(input_revel=args.input, input_bed=args.bed, output_tsv=args.output, output_ht=args.output_ht)
