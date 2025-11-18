import functools

import loguru

from cpg_utils import config
from metamist import graphql

LONG_READ_STRING = 'LongRead'
METAMIST_ANALYSIS_QUERY = graphql.gql(
    """
    query MyQuery($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
                output
                timestampCompleted
                meta
            }
        }
    }
""",
)


@functools.cache
def query_for_latest_analysis(
    dataset: str,
    stage_name: str | None = None,
) -> str | None:
    """
    Query for the latest analysis object of a given type in the requested project.

    Analysis entries for Talos all have unique types, so we can use this generic query method

    Args:
        dataset (str):         project to query for
        stage_name (str):      optional, if set, will only return entries with meta.stage == this
    Returns:
        str, the path to the latest object for the given type, or log a warning and return None
    """

    # quick escape here if I want to force usage of a specific MT (for testing, review later)
    if result := config.config_retrieve(['workflow', 'use_this_mt'], False):
        return result

    sequencing_type = config.config_retrieve(['workflow', 'sequencing_type'])
    long_read = config.config_retrieve(['workflow', 'long_read'], False)

    # swapping to a string we can freely modify
    query_dataset = dataset
    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    loguru.logger.info(f'Querying for MatrixTable in {query_dataset}')

    result = graphql.query(METAMIST_ANALYSIS_QUERY, variables={'dataset': query_dataset, 'type': 'matrixtable'})

    # get all the relevant entries, and bin by date
    analysis_by_date = {}
    for analysis in result['project']['analyses']:
        if analysis['output'] and (sequencing_type in {'all', analysis['meta'].get('sequencing_type')}):
            # skip over the partial-cohort AnnotateDataset objects
            if '_families-' in analysis['output']:
                loguru.logger.debug(
                    f'Skipping analysis {analysis["output"]} for dataset {query_dataset}. '
                    f'It is a partial-cohort AnnotateDataset object',
                )
                continue

            # manually implementing an XOR check - long read (bool) and LongRead in output must match
            if long_read != (LONG_READ_STRING in analysis['output']):
                loguru.logger.debug(
                    f'Skipping analysis {analysis["output"]} for dataset {query_dataset}. '
                    f'It does not match query parameter long_read={long_read}',
                )
                continue

            if stage_name is not None and analysis['meta'].get('stage') != stage_name:
                continue

            analysis_by_date[analysis['timestampCompleted']] = analysis['output']

    if not analysis_by_date:
        raise ValueError(f'No Analysis Entries found for dataset {query_dataset}')

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    return analysis_by_date[sorted(analysis_by_date)[-1]]
