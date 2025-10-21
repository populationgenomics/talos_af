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
            chrom = llist[0]
            start = int(llist[1])
            end = int(llist[2])
            bed_regions[chrom].append((start, end))

    return bed_regions


def region_of_interest(regions: REGION_DICT, chrom: str, pos: int) -> bool:
    """Check if this is a region of interest."""
    if chrom not in regions:
        return False

    for start, end in regions[chrom]:
        return start <= pos <= end
    return False
