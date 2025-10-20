#!/usr/bin/env python3

"""
Parses an ACMG specification and a GFF3 file. Generates a BED file of gene regions +/- padding, with the columns:
- chrom
- start
- end
"""

import gzip
import json
import re
import sys
from argparse import ArgumentParser

import loguru

CHROM_INDEX = 0
RESOURCE_INDEX = 1
TYPE_INDEX = 2
START_INDEX = 3
END_INDEX = 4
DETAILS_INDEX = 8

# +/- this is added to each gene region, this default can be overridden
FLANKING_REGION = 2000

# regular expressions to parse out sections of the GFF3 annotations
GENE_ID_RE = re.compile(r'gene:(ENSG\d+);')


def main(input_spec: str, gff3_file: str, output: str):
    """
    Read the GFF3 file, and generate a BED file of gene regions, plus padding

    Args:
        input_spec (str): the input GFF3 file
        gff3_file (str): path to the GFF3 file
        output (str): path to the intended BED output file
    """

    loguru.logger.info(f'Reading gff3 file from {gff3_file}, against specification {input_spec}')
    unmerged_lines = generate_bed_lines(specification_file=input_spec, gff3_file=gff3_file)
    merge_output(unmerged_lines=unmerged_lines, output=output)


def generate_bed_lines(
    specification_file: str,
    gff3_file: str,
) -> list[tuple[str, int, int]]:
    """
    Generate the new BED file, and return the lines as a list of lists for merging.
    """

    # read the file containing the specification for this round of analysis
    with open(specification_file, encoding='utf-8') as input_handle:
        specification = json.load(input_handle)

    output_lines: list[tuple[str, int, int]] = []

    # open and iterate over the GFF3 file
    with gzip.open(gff3_file, 'rt') as handle:
        for line in handle:
            # skip over headers and dividing lines
            if line.startswith('#'):
                continue

            line_as_list = line.rstrip().split('\t')

            # skip over non-genes (e.g. pseudogenes, ncRNA), only focus on Ensembl genes/transcripts
            if line_as_list[TYPE_INDEX] != 'gene' or 'ensembl' not in line_as_list[RESOURCE_INDEX]:
                continue

            # extract the gene ID from the details field
            if gene_id_match := GENE_ID_RE.search(line_as_list[DETAILS_INDEX]):
                gene_id = gene_id_match.group(1)
            else:
                loguru.logger.debug(f'Failed to extract gene name from {line_as_list[DETAILS_INDEX]}')
                continue

            # check if the current gene is one we're interested in
            if gene_id not in specification:
                loguru.logger.debug(f'Found a gene ID {gene_id} not in specification, skipping')
                continue

            # append the line to the output
            output_lines.append(
                (
                    f'chr{line_as_list[CHROM_INDEX]}',
                    int(line_as_list[START_INDEX]) - FLANKING_REGION,
                    int(line_as_list[END_INDEX]) + FLANKING_REGION,
                ),
            )
    return output_lines


def merge_output(
    unmerged_lines: list[tuple[str, int, int]],
    output: str,
):
    """
    Take each line, resolve overlapping regions, and write out to a new file.
    """
    contig = None
    start = None
    end = None
    with open(output, 'w') as handle:
        for this_chrom, this_start, this_end in unmerged_lines:
            if contig is None:
                contig = this_chrom
                start = this_start
                end = this_end
                continue

            # change in contig = write the previous block, and reset values
            if this_chrom != contig:
                handle.write(f'{contig}\t{start}\t{end}\n')
                contig = this_chrom
                start = this_start
                end = this_end
                continue

            # adjacent blocks close enough, merge them, but don't write
            if this_start - end < FLANKING_REGION:
                end = max(end, this_end)

            # far enough apart, write the previous block and reset
            else:
                handle.write(f'{contig}\t{start}\t{end}\n')
                start = this_start
                end = this_end

        # and write the final line
        handle.write(f'{contig}\t{start}\t{end}\n')


def cli_main():
    parser = ArgumentParser()
    parser.add_argument(
        '--input',
        help='Path to the ACMG criteria JSON specification',
        required=True,
    )
    parser.add_argument(
        '--gff3',
        help='Path to the compressed GFF3 file',
        required=True,
    )
    parser.add_argument(
        '--output',
        help='Path to output file, regions merged',
        required=True,
    )
    parser.add_argument('--debug', action='store_true', help='Debug logging mode')
    args = parser.parse_args()

    # opposite to the logging implementation, we have to strip the default and apply a non-debug log-level
    if not args.debug:
        loguru.logger.remove()
        loguru.logger.add(sys.stderr, level='INFO')

    main(
        input_spec=args.input,
        gff3_file=args.gff3,
        output=args.output,
    )


if __name__ == '__main__':
    cli_main()
