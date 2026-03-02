import gzip
from argparse import ArgumentParser
from typing import TYPE_CHECKING

import hail as hl
from loguru import logger

from cpg_utils import config, hail_batch, to_path

if TYPE_CHECKING:
    from hail.vds import VariantDataset


def get_reference_intervals() -> list[hl.Interval] | None:
    """
    Locates the bundled BED file and efficiently parses it into a list of Hail Intervals.
    """
    # 1. Locate the file
    interval_path = config.config_retrieve(['references', 'vds_read_intervals'], None)

    if interval_path is None:
        return None

    logger.info(f'Reading intervals from {interval_path}')

    intervals: list[hl.Interval] = []

    # 2. Parse it locally using standard Python
    with to_path(interval_path).open('rb') as f, gzip.open(f, 'rt') as gz:
        for line in gz:
            stripped = line.strip()
            # Skip empty and comment/header lines
            if not stripped or stripped.startswith(('#', 'track')):
                continue
            parts: list[str] = stripped.split('\t')
            if len(parts) < 3:
                logger.warning(f'Skipping malformed BED line: {line.rstrip()}')
                continue
            chrom, start, end = parts[:3]
            # Convert 0-based BED to 1-based string, add 1 to END as Hail is right exclusive by default
            start_locus = hl.Locus(chrom, int(start) + 1, reference_genome='GRCh38')
            end_locus = hl.Locus(chrom, int(end), reference_genome='GRCh38')
            # Convert 0-based BED to 1-based string, add 1 to END as Hail is right exclusive by default
            intervals.append(hl.Interval(start_locus, end_locus, includes_start=True, includes_end=True))

    return intervals


def main(input_vds: str, output_mt: str):
    """Load a sparse VariantDataset, export a dense MatrixTable with split multiallelics."""

    if qob_overrides := config.config_retrieve('qob_overrides'):
        hail_batch.init_batch(**qob_overrides)
    else:
        hail_batch.init_batch()

    vds: VariantDataset = hl.vds.read_vds(
        input_vds,
        intervals=get_reference_intervals(),
    )

    logger.info('Densifying data...')
    mt: hl.MatrixTable = hl.vds.to_dense_mt(vds)

    # remove gvcf_info if present
    if 'gvcf_info' in mt.row_value:
        mt = mt.drop('gvcf_info')

    logger.info('Post densification schema')
    mt.describe(handler=logger.info)

    # remove any monoallelic or non-alt-in-any-sample sites
    mt = mt.filter_rows((hl.len(mt.alleles) > 1) & (hl.agg.any(mt.LGT.is_non_ref())))

    # annotate with densified representations
    mt = mt.annotate_entries(
        GT=hl.vds.lgt_to_gt(mt.LGT, mt.LA),
        AD=hl.vds.local_to_global(
            mt.LAD,
            mt.LA,
            n_alleles=hl.len(mt.alleles),
            fill_value=0,
            number='R',
        ),
    )

    logger.info('Post Local-to-global reannotation schema')
    mt.describe(handler=logger.info)

    mt = mt.select_entries('GT', 'GQ', 'DP', 'AD')

    # split out multiallelic rows
    mt = hl.split_multi_hts(mt)

    logger.info('Post variant splitting schema')
    mt.describe(handler=logger.info)

    mt.write(output_mt, overwrite=True)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='Input VDS')
    parser.add_argument('--output', required=True, help='Output MatrixTable')
    args = parser.parse_args()
    main(input_vds=args.input, output_mt=args.output)
