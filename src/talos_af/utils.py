import json
import re
import sys
from collections import defaultdict
from functools import cache
from itertools import combinations_with_replacement
from typing import Any

import cyvcf2
import httpx
from cloudpathlib.anypath import to_anypath
from loguru import logger
from mendelbrot.pedigree_parser import PedigreeParser

from talos_af import models

REGION_DICT = dict[str, list[tuple[int, int]]]

HOMREF: int = 0
HETALT: int = 1
UNKNOWN: int = 2
HOMALT: int = 3

# in cyVCF2, these ints represent HOMREF, and UNKNOWN
BAD_GENOTYPES: set[int] = {HOMREF, UNKNOWN}
NON_HOM_CHROM = ['X', 'Y', 'MT', 'M']
X_CHROMOSOME = {'X'}

# used by cyvcf2 for... not phased
PHASE_SET_DEFAULT = -2147483648

GeneDict = dict[str, list[models.VariantAf]]
CompHetDict = dict[str, dict[str, list[models.VariantAf]]]
PATHOGENIC = 'Pathogenic/Likely Pathogenic'

CRITICAL_CSQ_DEFAULT = [
    'frameshift',
    'splice_acceptor',
    'splice_donor',
    'start_lost',
    'stop_gained',
    'stop_lost',
    'transcript_ablation',
]
CSQ_STRING = [
    'consequence',
    'gene',
    'transcript',
    'biotype',
    'strand',
    'amino_acid_change',
    'dna_change',
]
PHASE_BROKEN: bool = False

TRUNCATING = {'nonsense', 'frameshift'}

TYPE_RE = re.compile(r'p\.(?P<ref>\D)(?P<codon>\d+)(?P<alt>\D)')


def format_logger(
    log_level: int = logger.level('INFO').no,
    fmt_string: str = '{time:YYYY-MM-DD HH:mm:ss} - {file.path}:{line} - {level} - {message}',
    coloured: bool = False,
) -> None:
    """
    loguru is a cleaner interface than the standard logging module, but it doesn't allow for multiple instances
    instead of calling a get_logger function which returns a logger, we assume that any module using logging has
    imported `from loguru import logger` to get access to the logger.

    loguru.logger is also resistant to deepcopy, so there really is only a single global instance, meaning that the
    display/formatting of the logger is global to the entire process, and should only be set once.

    This helper method formats the logger instance with the given parameters, stripping out any previous handlers
    Because the global logger instance is modified, there is no return value

    >>> from loguru import logger
    >>> from cpg_flow.utils import format_logger
    >>> format_logger(log_level=10, fmt_string='{time} {level} {message}', coloured=True)
    >>> logger.info('This is an info message')

    Args:
        log_level (int): logging level, defaults to INFO. Can be overridden by config
        fmt_string (str): format string for this logger, defaults to DEFAULT_LOG_FORMAT
        coloured (bool): whether to colour the logger output
    """

    # Remove any previous loguru handlers
    logger.remove()

    # Add loguru handler with given format and level
    logger.add(
        sys.stdout,
        level=log_level,
        format=fmt_string,
        colorize=coloured,
        enqueue=True,
    )


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
            if not line.rstrip():
                continue
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

    return any(start <= pos <= end for start, end in regions[chrom])


def get_non_ref_samples(variant: 'cyvcf2.Variant', samples: list[str]) -> tuple[set[str], set[str]]:
    """
    For this variant, find all samples with a call. cyvcf2 uses 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT.
    Returns 2 sets of strings; het sample ID, hom sample IDs
    """
    het_samples = set()
    hom_samples = set()

    # this iteration is based on the cyvcf2 representations
    for sam, genotype_int in zip(samples, variant.gt_types, strict=True):
        if genotype_int in BAD_GENOTYPES:
            continue
        if genotype_int == HETALT:
            het_samples.add(sam)
        if genotype_int == HOMALT:
            hom_samples.add(sam)

    return het_samples, hom_samples


def get_phase_data(samples: list[str], var: 'cyvcf2.Variant') -> dict[str, dict[int, str]]:  # noqa: C901, PLR0912
    """
    read phase data from this variant

    we tolerate a variety of failures, the worst case scenario is that we have no phase data
    PS - PhaseSet
    PID/PGT - PhaseID/PhaseGenotype

    No other phase types are currently supported

    Args:
        samples (list[str]): all samples in the VCF
        var (cyvcf2.Variant):
    """

    global PHASE_BROKEN  # noqa: PLW0603

    phased_dict: dict[str, dict[int, str]] = defaultdict(dict)

    # check we have relevant attributes in the variant format fields
    if not any(x in var.FORMAT for x in ['PS', 'PGT', 'PID']):
        return dict(phased_dict)

    # first set the numpy.ndarray to be a list of ints
    # then zip against ordered sample IDs
    # this might need to store the exact genotype too
    # i.e. 0|1 and 1|0 can be in the same phase-set
    # but are un-phased variants
    try:
        if 'PS' in var.FORMAT:
            for sample, phase, genotype in zip(samples, map(int, var.format('PS')), var.genotypes, strict=True):
                # cyvcf2.Variant holds two ints, and a bool for biallelic calls
                # but only one int and a bool for hemi
                if len(genotype) == 3:
                    allele_1, allele_2, phased = genotype
                    gt = f'{allele_1}|{allele_2}'
                elif len(genotype) == 2:
                    allele_1, phased = genotype
                    gt = f'{allele_1}'
                else:
                    raise ValueError(f'Unexpected genotype length: {len(genotype)}: {genotype}')

                if not phased:
                    continue

                # phase set is a number
                if phase != PHASE_SET_DEFAULT:
                    phased_dict[sample][phase] = gt

        elif all(attr_name in var.FORMAT for attr_name in ['PGT', 'PID']):
            logger.info('Failed to find PS phase attributes')
            try:
                # retry using PGT & PID
                for sample, phase_gt, phase_id in zip(
                    samples,
                    var.format('PGT'),
                    var.format('PID'),
                    strict=True,
                ):
                    if phase_gt != '.' and phase_id != '.':
                        phased_dict[sample][phase_id] = phase_gt

            except KeyError as ke2:
                logger.info('Failed to determine phase information using PID and PGT')
                raise ke2
        elif not PHASE_BROKEN:
            logger.info('Found no PS phase attributes (known formats are PS, PGT/PID)')
            PHASE_BROKEN = True

    except (KeyError, ValueError):
        if not PHASE_BROKEN:
            logger.info('Failed to correctly parse known phase attributes using existing methods')
            logger.info(
                'Please post an issue on the Talos GitHub Repo with the VCF FORMAT lines and descriptions',
            )
        PHASE_BROKEN = True

    return dict(phased_dict)


def get_revel_as_dict(revel_details: str | None) -> dict[str, str]:
    if revel_details in {None, '.'}:
        return {}

    score, transcripts = revel_details.split('~')
    return dict.fromkeys(transcripts.split(','), score)


def organise_csq(
    var_details: dict[str, Any],
    symbol_indexed_lookup: dict[str, dict[str, str]],
) -> bool:
    """
    Read the transcript consequences, split into a list of details, integrate Revel and AlphaMissense where applicable.
    """

    bcsq = var_details.pop('bcsq', None)

    # no consequences?
    if not isinstance(bcsq, str):
        var_details['transcript_consequences'] = []
        return False

    # get and split revel data
    revel_dict = get_revel_as_dict(var_details.pop('revel', '.'))

    am_dict: dict[str, dict[str, str | float]] = {}
    if var_details.get('am_class', '.') != '.':
        transcript = str(var_details.pop('am_transcript'))
        am_class = str(var_details.pop('am_class'))
        am_score = float(var_details.pop('am_score'))
        am_dict[transcript] = {
            'am_class': am_class,
            'am_score': am_score,
        }
    else:
        var_details.pop('am_transcript', None)
        var_details.pop('am_score', None)

    consequential = False

    for txcsq in bcsq.split(','):
        # not zipping as Strict, sometimes the consequence is truncated
        txcsq_dict = dict(zip(CSQ_STRING, txcsq.split('|'), strict=False))

        if txcsq_dict['gene'] not in symbol_indexed_lookup:
            continue

        # skip over non-Mane transcripts
        if txcsq_dict['transcript'] != symbol_indexed_lookup[txcsq_dict['gene']]['enst']:
            continue

        txcsq_dict['ensg'] = symbol_indexed_lookup[txcsq_dict['gene']]['gene_id']

        if any(each_csq in CRITICAL_CSQ_DEFAULT for each_csq in txcsq_dict['consequence'].split('&')):
            consequential = True

        # integrate REVEL if appropriate
        if score := revel_dict.get(txcsq_dict['transcript']):
            txcsq_dict['revel'] = score

        # integrate alphamissense if appropriate
        if txcsq_dict['transcript'] in am_dict:
            txcsq_dict.update(**am_dict[txcsq_dict['transcript']])  # type: ignore[arg-type]

        var_details['transcript_consequences'] = txcsq_dict

    return consequential


def create_small_variant(
    var: 'cyvcf2.Variant',
    samples: list[str],
    symbol_indexed_lookup: dict[str, dict[str, str]],
) -> models.VariantAf | None:
    """Takes a small variant and creates a Model from it."""

    coordinates = models.Coordinates(chrom=var.CHROM.replace('chr', ''), pos=var.POS, ref=var.REF, alt=var.ALT[0])

    # this is so flexible it's basically useless
    info: dict[str, Any] = {x.lower(): y for x, y in var.INFO if y is not None}

    clinvar_path = info.get('clinical_significance') == PATHOGENIC

    # we never want the info dict to contain None values, so only insert clinvar entries if pathogenic
    if clinvar_path:
        info |= {
            'clinvar_allele': info['allele_id'],
            'clinvar_stars': int(info['gold_stars']),
        }

    consequential = organise_csq(info, symbol_indexed_lookup)

    if not (clinvar_path or consequential):
        return None

    het_samples, hom_samples = get_non_ref_samples(variant=var, samples=samples)

    phased = get_phase_data(samples, var)

    # TODO resolve this at some point so we can have non-transcript ClinVar Path vars come through
    if 'transcript_consequences' not in info:
        return None

    transcript_consequences: dict[str, str | int | float] = info.pop('transcript_consequences')

    # revise in the future, for now we're skipping any non-transcript consequence variants
    if not transcript_consequences:
        return None

    return models.VariantAf(
        gene=transcript_consequences['ensg'],
        coordinates=coordinates,
        clinvar_path=clinvar_path,
        high_impact=consequential,
        info=info,
        het_samples=het_samples,
        hom_samples=hom_samples,
        phased=phased,
        transcript_consequences=transcript_consequences,
    )


def gather_gene_dict_from_vcf(
    vcf_path: str,
    acmg_spec: dict[str, dict[str, str]],
) -> dict[str, list[models.VariantAf]]:
    """
    takes a cyvcf2.VCFReader instance, and a specified chromosome
    iterates over all variants in the region, and builds a lookup

    optionally takes a second VCF and incorporates into same dict

    Args:
        vcf_path (str): the VCF to read
        acmg_spec (dict): the ACMG specification

    Returns:
        A lookup in the form
        {
            gene1: [var1, var2],
            gene2: [var3],
            ...
        }
    """

    # a dict to allow lookup of variants on this whole chromosome
    variant_count = 0
    skipped_variant_count = 0
    contig_dict = defaultdict(list)

    variant_source = cyvcf2.VCF(vcf_path)

    symbol_indexed_acmg = {value['gene']: value for key, value in acmg_spec.items()}

    # iterate over all variants on this contig and store by unique key
    # if contig has no variants, prints an error and returns []
    for variant_row in variant_source:
        variant = create_small_variant(
            var=variant_row, samples=variant_source.samples, symbol_indexed_lookup=symbol_indexed_acmg
        )
        if variant is None:
            skipped_variant_count += 1
            continue

        # update the variant count
        variant_count += 1
        # update the gene index dictionary
        contig_dict[variant.gene].append(variant)

    logger.info(
        f"""
        VCF {vcf_path}:
        \t{variant_count} selected variants
        \t{skipped_variant_count} skipped variants
        \tin {len(contig_dict)} genes
        """
    )

    return dict(contig_dict)


def find_comp_hets(var_list: list[models.VariantAf], pedigree: PedigreeParser) -> CompHetDict:
    """
    Find compound het pairs, variants provided in the format [var1, var2, ...]

    generates pair content in the form
    {
        sample: {
            var_as_string: [partner_variant, ...],
        }
    }
    """

    # create an empty dictionary
    comp_het_results: CompHetDict = defaultdict(dict)

    # use combinations_with_replacement to find all gene pairs
    for var_1, var_2 in combinations_with_replacement(var_list, 2):
        if (var_1.coordinates == var_2.coordinates) or var_1.coordinates.chrom in NON_HOM_CHROM:
            continue

        # iterate over any samples with a het overlap
        for sample in var_1.het_samples.intersection(var_2.het_samples):
            phased = False

            # don't assess male compound hets on sex chromosomes
            if pedigree.participants[sample].sex == 1 and var_1.coordinates.chrom in X_CHROMOSOME:
                continue

            # check for both variants being in the same phase set
            if sample in var_1.phased and sample in var_2.phased:
                # check for presence of the same phase set
                for phase_set in [ps for ps in var_1.phased[sample] if ps in var_2.phased[sample]]:
                    if var_1.phased[sample][phase_set] == var_2.phased[sample][phase_set]:
                        phased = True
            if not phased:
                comp_het_results[sample].setdefault(var_1.coordinates.string_format, []).append(var_2)
                comp_het_results[sample].setdefault(var_2.coordinates.string_format, []).append(var_1)

    return comp_het_results


def is_variant_truncating(variant: models.VariantAf) -> bool:
    """Check if any of the transcript consequences for this variant are truncating."""
    for each_txc in variant.transcript_consequences:
        # split out separate terms so we can compare separately
        csq_set = each_txc['consequence'].split('&')
        if any(truncating_term in csq_set for truncating_term in TRUNCATING):
            return True
    return False


def is_variant_exact_p(change: str, variant: models.VariantAf) -> bool:
    """
    Check if any of the phased transcript consequences for this variant are exact protein changes.
    Because I'm using SamTools for the annotation I've got to revert the schema changes to that AA change notation.
    """
    match = re.match(TYPE_RE, change)
    codon = match.group('codon')
    ref = match.group('ref')
    alt = match.group('alt')
    samtools_format_string = f'{codon}{ref}>{codon}{alt}'

    return any(each_txc['amino_acid_change'] == samtools_format_string for each_txc in variant.transcript_consequences)


def apply_gene_specific_rules(
    rule: str,
    variant: models.VariantAf,
) -> bool:
    """
    Take the initial list of variants we selected, and apply more specific rules to back-filter.
    This is a bit janky and manual for now. To be reviewed.
    """
    if rule == 'truncating':
        return is_variant_truncating(variant)

    if rule.startswith('p.'):
        return is_variant_exact_p(change=rule, variant=variant)

    raise NotImplementedError('Not sure what the rule should be here.')


@cache
def get_latest_zenodo_file(record_id: int) -> tuple[str, str]:
    """Take a zenodo record ID, find the latest version of that record, and return the attachment link."""

    manual_retries = 3

    # use this record ID to find the latest record ID
    latest_url = f'https://zenodo.org/api/records/{record_id}/versions/latest'

    data: dict | None = None
    while manual_retries > 0:
        manual_retries -= 1
        try:
            response = httpx.get(latest_url, headers={'Accept': 'application/json'}, timeout=60, follow_redirects=True)
            data = response.json()
        except (httpx.ConnectTimeout, httpx.DecodingError):
            logger.error(f'Failed to get latest zenodo record ID: {record_id}')
            continue

    if data is None:
        logger.error(f'Failed to get latest zenodo record ID: {record_id}, exhausted retries')

    # navigate to the files, and find the content URL we can use to download
    file_link = data['files'][0]['links']['self']
    return str(data['recid']), file_link


def update_dates_from_prior_data(current: models.ResultsAf, prior_path: str | None = None) -> None:
    """
    Ingest a previous object, and update any current discovery dates if they were previously seen

    Args:
        current: models.ResultsAf object, containing this run's results
        prior_path: str, path to prior data
    """
    if prior_path is None:
        return

    with to_anypath(prior_path).open() as handle:
        json_data = json.load(handle)
        # potentially walk-up model version
        old_data = models.ResultsAf.model_validate(json_data)

        # iterate over the older object
        for sample, old_instances in old_data.instances.items():
            # if the older data sample ID is not in the current dataset, move on
            if sample not in current.instances:
                continue

            # find the associated date for all the variant instances for this sample ID
            old_id_dates = {instance.var_id: instance.first_seen for instance in old_instances[sample]}
            for new_instance in current.instances[sample]:
                # if the same var_id is seen in this current round, backdate the first_seen date
                if new_instance.var_id in old_id_dates:
                    new_instance.first_seen = old_id_dates[new_instance.var_id]
