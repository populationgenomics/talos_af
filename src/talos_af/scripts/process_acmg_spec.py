"""
script for taking the ACMG specification

This is aimed at the supplementary table from the publication:
 - https://www.sciencedirect.com/science/article/pii/S1098360025001017#mmc1
 - this is natively in a xlsx format, which is harder to parse
 - open the spreadsheet, trim off the top two rows, and remove any trailing non-data rows
 - save the result as a TSV, and use that as input to this script

The results of running this script on criteria v3.3 will be enclosed here.

The other required file is MANE.GRCh38.v1.4.summary.txt.gz from the NCBI FTP site:
 - https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz
 - the ACMG release only contains gene symbols, so we use this to map these to ENSG IDs for each gene
 - the MANE file has the genomic regions for each gene, albeit with NC_nnn.n contig numbering
"""

import gzip
import json
import re
from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader
from functools import cache

import loguru

# ordered list of contigs
CONTIG_SORT_ORDER = [f'chr{x}' for x in [*list(range(1, 23)), 'X', 'Y']]
CONTIG_START = re.compile(r'NC_0+')

MANE_PLUS_CLINICAL = 'MANE Plus Clinical'
MANE_SELECT = 'MANE Select'
P_AND_LP = 'P and LP'
SPECIFIC_CHANGE = re.compile(r'(p.[A-Z][0-9]+[A-Z]) (?:\w+ )?only')

# define padding params
GENE_PAD = 500  # margin to add at start/end for each gene
BRIDGING_PAD = 2000  # if genes are within this span, merge regions


@cache
def nc_to_chr(nc: str) -> str:
    """Takes a NC contig string, and converts it to a chrN format contig. Must be a better solution to this."""

    # strip off version
    nc_start = nc.split('.')[0]

    # strip off leading digits
    leading_digits = re.match(CONTIG_START, nc_start).group()
    nc_tail = nc_start.replace(leading_digits, '')

    # now return the nth contig
    return CONTIG_SORT_ORDER[int(nc_tail) - 1]


def get_mane_mapping(mane_path: str) -> dict:
    """
    Relevant column headings:
    - Ensembl_Gene: ENSG ID with version
    - symbol: Gene symbol
    - RefSeq_nuc: NM ID with version
    - Ensembl_nuc: ENST with version
    - GRCh38_chr: NC_nnn.n contig numbering
    - chr_start: int
    - chr_end: int

    Returns:
        Dict, indexed on symbol, containing a corresponding ENSG, and a preferred refseq & ensembl transcript ID
    """

    mane_mapping: dict = {}
    with gzip.open(mane_path, 'rt') as mane_file:
        mane_reader = DictReader(mane_file, delimiter='\t')
        for line in mane_reader:
            symbol = line['symbol']
            ensg = line['Ensembl_Gene'].split('.')[0]
            nm_id = line['RefSeq_nuc'].split('.')[0]
            enst = line['Ensembl_nuc'].split('.')[0]
            ensp = line['RefSeq_prot'].split('.')[0]

            # allow for MANE Select and Plus Clinical to be populated against the same gene ID/symbol

            if not line['GRCh38_chr'].startswith('NC'):
                loguru.logger.debug(f'Skipping {symbol} as non-canonical contig ({line["GRCh38_chr"]})')
                continue

            contig = nc_to_chr(line['GRCh38_chr'])

            # if we found a previous row, add in a new MANE type entry
            if symbol in mane_mapping:
                mane_mapping[symbol]['mane'][line['MANE_status']] = {
                    'nm_id': nm_id,
                    'enst': enst,
                    'ensp': ensp,
                }

            # first row, add the whole record
            else:
                mane_mapping[symbol] = {
                    'ensg': ensg,
                    'mane': {
                        line['MANE_status']: {
                            'nm_id': nm_id,
                            'enst': enst,
                            'ensp': ensp,
                        },
                    },
                    'contig': contig,
                    'start': line['chr_start'],
                    'end': line['chr_end'],
                }
    return mane_mapping


def eat_acmg_tsv(acmg_file: str) -> dict[str, dict[str, str]]:
    """Smash that ACMG spec file."""
    acmg_digest: dict[str, dict[str, str]] = defaultdict(dict)

    with open(acmg_file) as input_file:
        reader = DictReader(input_file, delimiter='\t')
        for row in reader:
            gene = row['Gene'].rstrip()

            row_dict = {
                'gene': gene,
                'moi': row['Inheritance'].rstrip(),
            }

            vars_to_report = row['Variants to report']

            if 'truncating' in vars_to_report:
                row_dict['reportable'] = 'consequence'
                row_dict['specific_type'] = 'truncating'

            elif P_AND_LP in vars_to_report:
                row_dict['reportable'] = 'all'

            elif result := re.findall(SPECIFIC_CHANGE, vars_to_report):
                row_dict['reportable'] = 'specific'
                row_dict['specific_type'] = result[0]

            acmg_digest[gene] = row_dict

    return acmg_digest


def get_merged_bed_regions(acmg_digest: dict[str, dict[str, str]]) -> list[tuple[str, int, int]]:
    """Pick out all the locations, add some padding, sort and merge. EZ."""

    all_regions: list[tuple[str, int, int]] = []

    # pick up all the regions, indexed by contig (ACMG default ordering is alphabetical by gene symbol)
    contig_regions: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for acmg_section in acmg_digest.values():
        contig_regions[acmg_section['contig']].append(
            (int(acmg_section['start']) - GENE_PAD, int(acmg_section['end']) + GENE_PAD)
        )

    for contig in CONTIG_SORT_ORDER:
        if contig not in contig_regions:
            continue

        last_start: int | None = None
        last_end: int | None = None

        # get the sorted list of regions, ordered by region start
        regions = sorted(contig_regions[contig], key=lambda couple: couple[0])
        for region in regions:
            if last_end is None:
                last_start = region[0]
                last_end = region[1]
                continue
            # if the gap between this and previous region is small, merge them
            if last_end + BRIDGING_PAD > region[0]:
                last_end = max(last_end, region[1])
            # otherwise add the region to the list and grab new values
            else:
                all_regions.append((contig, last_start, last_end))
                last_start = region[0]
                last_end = region[1]
        # end of contig, drop it onto the list
        all_regions.append((contig, last_start, last_end))

    return all_regions


def main(input_spec: str, json_out: str, bed_out: str, mane_file: str) -> None:
    acmg_digest = eat_acmg_tsv(input_spec)

    mane_mapping = get_mane_mapping(mane_file)

    # update the acmg dict to contain everything
    for acmg_key, acmg_values in acmg_digest.items():
        acmg_values |= mane_mapping[acmg_key]  # noqa: PLW2901

    all_regions = get_merged_bed_regions(acmg_digest)

    with open(bed_out, 'w') as handle:
        for contig, start, end in all_regions:
            handle.write(f'{contig}\t{start}\t{end}\n')

    with open(json_out, 'w', encoding='utf-8') as output_file:
        json.dump(acmg_digest, output_file, indent=4)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='ACMG specification csv file')
    parser.add_argument('--mane', help='The MANE text file')
    parser.add_argument('--json_out', help='Path to write a minimised spec to.')
    parser.add_argument('--bed_out', help='Path to write a bed file out to.')
    args = parser.parse_args()
    main(input_spec=args.input, bed_out=args.bed_out, json_out=args.json_out, mane_file=args.mane)
