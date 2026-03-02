from argparse import ArgumentParser
from typing import TYPE_CHECKING

import hail as hl
from loguru import logger

from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hail.vds import VariantDataset


def get_reference_intervals() -> list[hl.Interval] | None:
    """
    Locates the bundled BED file and efficiently parses it into a list of Hail Intervals.
    """
    if interval_path := config.config_retrieve(['references', 'vds_read_intervals'], None):
        logger.info(f'Reading intervals from {interval_path}')
        interval_table = hl.import_bed(interval_path, force_bgz=True)
        return interval_table.interval.collect()
    return None


def main(input_vds: str, output_mt: str):
    """Load a sparse VariantDataset, export a dense MatrixTable with split multiallelics."""

    hail_batch.init_batch()

    intervals_for_reading = get_reference_intervals()

    vds: VariantDataset = hl.vds.read_vds(input_vds, intervals=intervals_for_reading)

    logger.info('Densifying data...')
    mt: hl.MatrixTable = hl.vds.to_dense_mt(vds)

    # taken from _filter_rows_and_add_tags in large_cohort/site_only_vcf.py
    # remove any monoallelic or non-ref-in-any-sample sites
    mt = mt.filter_rows((hl.len(mt.alleles) > 1) & (hl.agg.any(mt.LGT.is_non_ref())))

    mt = hl.experimental.sparse_split_multi(mt)

    mt = mt.select_entries('GT', 'GQ', 'DP', 'AD')

    mt.write(output_mt, overwrite=True)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='Input VDS')
    parser.add_argument('--output', required=True, help='Output MatrixTable')
    args = parser.parse_args()
    main(input_vds=args.input, output_mt=args.output)
