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
"""

import gzip
import json
import re
from argparse import ArgumentParser
from csv import DictReader

P_AND_LP = 'P and LP'
SPECIFIC_CHANGE = re.compile(r'(p.[A-Z][0-9]+[A-Z]) (?:\w+ )?only')


def get_mane_mapping(mane_path: str) -> dict[str, dict[str, str]]:
    """
    Relevant column headings:
    - Ensembl_Gene: ENSG ID (with trailing decimal)
    - symbol: Gene symbol
    - RefSeq_nuc: NM ID with trailing decimal
    - Ensembl_nuc: ENST with trailing decimal

    Returns:
        Dict, indexed on symbol, containing a corresponding ENSG, and a preferred refseq & ensembl transcript ID
    """

    mane_mapping = {}
    with gzip.open(mane_path, 'rt') as mane_file:
        mane_reader = DictReader(mane_file, delimiter='\t')
        for line in mane_reader:
            symbol = line['symbol']
            ensg = line['Ensembl_Gene'].split('.')[0]
            nm_id = line['RefSeq_nuc'].split('.')[0]
            enst = line['Ensembl_nuc'].split('.')[0]
            mane_mapping[symbol] = {
                'ensg': ensg,
                'nm_id': nm_id,
                'enst': enst,
            }
    return mane_mapping


def main(input_spec: str, output_path: str, mane_file: str) -> None:
    parsed: dict = {}
    mane_mapping = get_mane_mapping(mane_file)

    with open(input_spec) as input_file:
        reader = DictReader(input_file, delimiter='\t')
        for row in reader:
            gene = row['Gene']

            # hard lookup here, we don't expect or allow failures. These are super well known genes
            mane_content = mane_mapping[gene]

            row_dict = {
                'gene': gene,
                'moi': row['Inheritance'],
                'gene_id': mane_content['ensg'],
                'nm_id': mane_content['nm_id'],
                'enst': mane_content['enst'],
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

            parsed[mane_content['ensg']] = row_dict

    with open(output_path, 'w', encoding='utf-8') as output_file:
        json.dump(parsed, output_file, indent=4)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', help='ACMG specification csv file')
    parser.add_argument('--output', help='Location to write a minimised spec to')
    parser.add_argument('--mane', help='The MANE text file')
    args = parser.parse_args()
    main(args.input, args.output, args.mane)
