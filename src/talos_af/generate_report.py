"""
Methods for taking the final output and generating static report content

This is another total rewrite, which tries to fit some resource-friendly
frontage onto the report, so that it loads in good time.

If there's a common prefix (e.g. by year), we split the data into sub-reports,
but we don't need to keep paring it down and down
"""

import json
import re
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import jinja2
import pandas as pd
from cloudpathlib.anypath import to_anypath
from loguru import logger
from mendelbrot.bcftools_interpreter import TYPES_RE, classify_change

from talos_af.config import config_retrieve
from talos_af.models import ReportableVariant, ResultsAf
from talos_af.utils import format_logger

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'

# above this length we trim the actual bases to just an int
MAX_INDEL_LEN: int = 10

# regex pattern - number, number, not-number
CDNA_SQUASH = re.compile(r'(?P<type>ins|del)(?P<bases>[ACGT]+)$')
MEAN_SLASH_SAMPLE = 'Mean/sample'


def parse_ids_from_file(ext_id_file: str | None) -> dict[str, str] | None:
    """
    Reads a file containing external IDs and returns a dictionary mapping
    this can be headerless TSV, CSV, in which case the two columns are Sample ID (in VCF/MT) and an external ID
    If provided in JSON format, this must be a dictionary of Sample ID -> External ID
    """

    # escape if there was nothing provided, or the file doesn't exist
    if ext_id_file is None:
        return None
    if not to_anypath(ext_id_file).exists():
        logger.warning(f'External ID file {ext_id_file} does not exist or was not accessible, skipping')
        return None

    id_mapping: dict[str, str] = {}

    file_as_path = to_anypath(ext_id_file)
    if (suffix := file_as_path.suffix) == '.json':
        # read the JSON file as a dictionary
        with file_as_path.open() as f:
            id_mapping = json.load(f)

        if not isinstance(id_mapping, dict):
            raise ValueError(f'Expected a dictionary in {ext_id_file}, got {type(id_mapping)}')

    elif suffix in ['.tsv', '.csv']:
        delimiter = ',' if suffix == '.csv' else '\t'
        with file_as_path.open(encoding='utf-8') as handle:
            for line in handle:
                # skip empty lines
                if not line.strip():
                    continue
                # split the line into two parts
                parts = line.strip().split(delimiter)
                if len(parts) != 2:
                    raise ValueError(f'Expected two columns in {ext_id_file}, got {len(parts)}')
                sample_id, ext_id = parts
                id_mapping[sample_id] = ext_id

    else:
        raise ValueError(f'Unsupported file format: {ext_id_file} - expected JSON, TSV or CSV, with correct extension.')

    return id_mapping


class NoVariantsFoundError(Exception):
    """raise if a report subset contains no data"""


@dataclass
class DataTable:
    """
    Representation of a DataTables table that the Jinja2 templating system renders.
    """

    id: str
    columns: list[str]
    rows: list[Any]
    heading: str = ''
    description: str = ''


class HTMLBuilder:
    """
    Takes the input, makes the output
    """

    def __init__(
        self,
        results_object: ResultsAf,
        link_engine: 'LinkEngine | None' = None,
        ext_id_map: dict[str, str] | None = None,
        igv_dir: str | None = None,
    ):
        """
        Args:
            results_object (ResultsAf): the results object
            link_engine (LinkEngine, optional): the link engine to generate hyperlinks with
            ext_id_map (dict[str, str], optional): a mapping of sample IDs to external IDs, optional
            igv_dir (str, optional): the igv_dir name to use to find CRAM, optional
        """

        if ext_id_map is None:
            logger.info('No External IDs were provided, using Sample names as IDs')

        self.ext_id_map = ext_id_map or {}

        # take the link-generating instance (can be None)
        self.link_engine = link_engine

        self.metadata = results_object.metadata

        # Process samples and variants
        self.samples: list[Sample] = []

        for sample, variant_instances in results_object.instances.items():
            if not variant_instances:
                continue

            sample_variants = [
                Variant(
                    # the annotations
                    results_object,
                    # this individual event
                    instance,
                )
                for instance in variant_instances
            ]
            self.samples.append(
                Sample(
                    name=sample,
                    variants=sample_variants,
                    igv_dir=igv_dir,
                ),
            )
        self.samples.sort(key=lambda x: x.name)

    @staticmethod
    def mane_to_string(mane_data: dict[str, dict[str, str]]) -> str:
        """Munge the MANE data block to be suitable for HTML embedding."""
        string_bits: list[str] = []
        for mane_type, mane_entry in mane_data.items():
            details = ', '.join(value for key, value in mane_entry.items())
            string_bits.append(f'{mane_type}: {details}')

        return '. '.join(string_bits)

    def read_metadata(self) -> dict[str, pd.DataFrame]:
        """
        parses into a general table
        """

        return {
            'Meta': pd.DataFrame(
                {
                    'Data': key.capitalize(),
                    'Value': self.metadata.__getattribute__(key),
                }
                for key in ['run_date']
            ),
            'Specification': pd.DataFrame(
                {
                    'Symbol': section.gene,
                    'MOI': section.moi,
                    'MANE': self.mane_to_string(section.mane),
                    'Reportable': section.reportable,
                    'Specific Type': section.specific_type,
                }
                for section in self.metadata.specification.values()
            ),
        }

    def write_html(self, output_filepath: str):
        """
        Uses the results to create the HTML tables
        writes all content to the output path

        Args:
            output_filepath (str): where to write the results to
        """

        meta_tables_raw = self.read_metadata()
        meta_tables = {
            name: {
                'columns': table.columns.tolist(),
                'rows': table.to_dict(orient='records'),
            }
            for name, table in meta_tables_raw.items()
            if not table.empty
        }

        # ensure Meta (run metadata) appears first and has a descriptive heading
        if 'Meta' in meta_tables:
            meta_tables = {'Run Metadata': meta_tables.pop('Meta')} | meta_tables

        template_context = {
            'run_date': self.metadata.run_date,
            'samples': self.samples,
            'report_title': 'Talos AF Report',
            'meta_tables': meta_tables,
            'config_json': json.dumps(config_retrieve([]), indent=2, sort_keys=True),
        }

        # write all HTML content to the output file in one go
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR),
            autoescape=True,
            trim_blocks=True,
            keep_trailing_newline=False,
            lstrip_blocks=True,
        )
        template = env.get_template('index.html.jinja')
        content = template.render(**template_context)
        with open(output_filepath, 'w') as handle:
            handle.writelines(content)

        logger.info(f'Wrote {output_filepath}')


class Sample:
    """
    Sample related logic
    """

    def __init__(
        self,
        name: str,
        variants: list['Variant'],
        igv_dir: str | None = None,
    ):
        self.name = name
        self.ext_id = 'EXTERNAL'
        self.variants = variants
        self.igv_dir = igv_dir

    def __str__(self):
        return self.name


class LinkEngine:
    """
    Generate links to external resources based on configuration settings
    """

    def __init__(
        self,
        template: str,
        variant_template: str | None = None,
        external: bool = False,
        lookup: str | None = None,
    ):
        """

        Args:
            template (str): mandatory - without this there's no sense generating an instance
            variant_template (str): optional, if there's a string here, we'll try and generate variant-specific links
            external (bool): if True, embed/lookup external IDs in the lookup dictionary. Default is sample.name.
                             This is mostly for a CPG internal use-case, where the seqr lookup and external lookup come
                             from different sources. The Lookup variable makes this redundant.
            lookup (dict): optional, a path to a CSV/TSV/JSON file, used to connect sample ID -> arbitrary ID
        """
        self.template = template
        self.variant_template = variant_template
        self.external = external
        self.lookup = parse_ids_from_file(lookup)

    def get_string_id(self, sample: Sample) -> str | None:
        """Get the string ID for the sample to use in links."""

        key = sample.ext_id if self.external else sample.name

        if self.lookup:
            # bail here instead of generating broken links
            if key not in self.lookup:
                return None

            return self.lookup[key]

        return key

    def generate_sample_link(self, sample: Sample):
        """Generates a sample/family level link using the template."""

        string_id = self.get_string_id(sample)

        # escape here - if we want an ID translated, don't return a hyperlink
        # feels better than returning a broken hyperlink
        if string_id is None:
            return None

        return self.template.format(sample=string_id)

    def generate_variant_link(self, sample: Sample, var_string: str) -> str | None:
        """Generate a Sample & Variant level link using the template."""

        if not self.variant_template:
            return None

        string_id = self.get_string_id(sample)

        # escape here - if we want an ID translated, don't return a hyperlink
        # feels better than returning a broken hyperlink
        if string_id is None:
            return None

        return self.variant_template.format(sample=self.get_string_id(sample), variant=var_string)


class Variant:
    """
    Handle as much of per variant logic as we can here. Hopefully, this is just simple
    data munging and mapping operations.

    Try not to put presentation logic here - keep it in the jinja templates
    """

    def get_var_change(self) -> str:
        """
        Find the variant change for the variant
        - we want to truncate huge small variant InDels (ballooning column width)
           - e.g. LOLOLOLOLOLOLOLOLOLOLOLOLOLOLOLOLO->A -> del 34bp
        - SVs always need to be presented differently
           - e.g INS 4079bp
        """
        if len(self.ref) > MAX_INDEL_LEN or len(self.alt) > MAX_INDEL_LEN:
            ref_len = len(self.ref)
            alt_len = len(self.alt)
            if ref_len > alt_len:
                return f'del {ref_len - alt_len}bp'
            if ref_len == alt_len:
                return f'complex delins {ref_len}bp'
            return f'ins {alt_len - ref_len}bp'

        return f'{self.ref}->{self.alt}'

    def __init__(self, report_object: ResultsAf, instance: ReportableVariant):
        # pick out the annotations for this one variant
        self.var_id = instance.var_id

        variant_annotations = report_object.variants[self.var_id]

        self.gene: str = variant_annotations.gene

        self.gt = instance.genotype

        self.acmg = report_object.metadata.specification[variant_annotations.gene]

        self.symbol = self.acmg.gene

        self.high_impact = variant_annotations.high_impact
        self.clinvar_path = variant_annotations.clinvar_path

        self.info = variant_annotations.info
        self.coordinates = variant_annotations.coordinates
        self.chrom = variant_annotations.coordinates.chrom
        self.pos = variant_annotations.coordinates.pos
        self.ref = variant_annotations.coordinates.ref
        self.alt = variant_annotations.coordinates.alt
        self.change = self.get_var_change()

        self.first_tagged: str = instance.first_seen
        self.support_vars = ', '.join(sorted(instance.support_vars))
        self.transcript_consequences = variant_annotations.transcript_consequences

        # Summarise CSQ strings per MANE entity
        self.mane_csq = self.parse_csq()

        # shove in the AM and REVEL scores, even if they're missing
        all_am_scores = [
            x.get('am_score') for x in self.transcript_consequences if isinstance(x.get('am_score'), float)
        ]
        all_revel_scores = [x.get('revel') for x in self.transcript_consequences if isinstance(x.get('revel'), float)]
        self.info |= {
            'am_score': max(all_am_scores) if all_am_scores else 'missing',
            'revel': max(all_revel_scores) if all_revel_scores else 'missing',
        }

    def __str__(self) -> str:
        return self.var_id

    def parse_csq(self) -> list[tuple[str, str, str]]:
        """
        Parse CSQ variant string returning:
            - set of "consequences" from MANE transcripts
            - Set of variant effects in p. nomenclature (or c. if no p. is available)

        condense massive cdna annotations, e.g.
        c.4978-2_4978-1insAGGTAAGCTTAGAAATGAGAAAAGACATGCACTTTTCATGTTAATGAAGTGATCTGGCTTCTCTTTCTA
        """

        mane_results: list[tuple[str, str, str]] = []

        # allow for a separate set of consequences per MANE Select/Plus Clinical
        for consequence_dict in self.transcript_consequences:
            consequences = set()
            nmd_consequences = set()

            consequence = consequence_dict.get('consequence')

            # for the types I know how to parse, update them
            if (type_match := TYPES_RE.match(consequence)) and consequence_dict.get('amino_acid_change'):
                consequence_dict['amino_acid_change'] = classify_change(
                    consequence_dict['amino_acid_change'],
                    consequence=type_match[0],
                )

            p_change = (
                f'{consequence_dict["transcript"]} - {consequence_dict["amino_acid_change"]}'
                if consequence_dict.get('amino_acid_change')
                else ''
            )

            csq_replaced = consequence.replace('_variant', '').replace('_', ' ')
            variant_csqs = csq_replaced.split('&')

            if 'NMD transcript' in variant_csqs:
                variant_csqs.remove('NMD transcript')
                nmd_consequences.update(variant_csqs)
            else:
                consequences.update(variant_csqs)

            # simplify the consequence strings
            if consequences:
                csq_string = ', '.join(consequences)
            elif nmd_consequences:
                csq_string = 'NMD only: ' + ', '.join(nmd_consequences)
            else:
                csq_string = 'None?'

            mane_results.append(
                (
                    consequence_dict['mane_type'],
                    csq_string,
                    p_change,
                ),
            )

        return mane_results


def cli_main():
    logger.info('Running HTML builder')
    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to analysis results', required=True)
    parser.add_argument('--output', help='Final HTML filename', required=True)
    parser.add_argument('--igv_dir', help='Used to generate IGV desktop links', default=None)
    parser.add_argument('--ext_ids', help='Optional, Mapping file for external IDs', default=None)
    parser.add_argument('--seqr_ids', help='Optional, Mapping file for Seqr IDs', default=None)
    args = parser.parse_args()
    main(
        results=args.input,
        output=args.output,
        igv_dir=args.igv_dir,
        ext_id_file=args.ext_ids,
        seqr_id_file=args.seqr_ids,
    )


def main(
    results: str,
    output: str,
    igv_dir: str,
    ext_id_file: str | None = None,
    seqr_id_file: str | None = None,
):
    """

    Args:
        results (str): path to the MOI-tested results file
        output (str): where to write the HTML file
        igv_dir (str): CPG-internal, used to identify location of a sample CRAM file
        ext_id_file (str | None): optional, path to a file containing external IDs
        seqr_id_file (str | None): optional, path to a file containing Seqr IDs
    """

    with to_anypath(results).open(mode='r') as f:
        f_json = json.load(f)
        results_object = ResultsAf.model_validate(f_json)

    # can be None if absent, or is a lookup of sample ID in VCF ~ an external ID
    external_id_map = parse_ids_from_file(ext_id_file)

    # set up the link builder, or None
    if (link_section := config_retrieve(['CreateTalosHTML', 'hyperlinks'], None)) and seqr_id_file:
        link_builder = LinkEngine(**link_section, lookup=seqr_id_file)
    elif link_section:
        link_builder = LinkEngine(**link_section, lookup=None)
    else:
        link_builder = None

    html = HTMLBuilder(
        results_object=results_object,
        link_engine=link_builder,
        ext_id_map=external_id_map,
        igv_dir=igv_dir,
    )

    # if this fails with a NoVariantsFoundException, there were no variants to present in the whole cohort
    # catch this, but fail gracefully so that the process overall is a success
    try:
        logger.debug(f'Writing whole-cohort categorised variants to {output}')
        html.write_html(output_filepath=output)
    except NoVariantsFoundError:
        logger.warning('No Categorised variants found in this whole cohort')


if __name__ == '__main__':
    format_logger(log_level=logger.level('INFO').no)
    cli_main()
