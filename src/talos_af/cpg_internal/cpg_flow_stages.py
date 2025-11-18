from argparse import ArgumentParser

from cpg_flow import stage, targets, workflow
from cpg_utils import Path, config, hail_batch

from talos_af.cpg_internal import utils as internal_utils


@stage.stage
class ParseAcmgSpec(stage.MultiCohortStage):
    """
    Read the TSV representing the ACMG specification, and parse it
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {'json': self.prefix / config.config_retrieve(['acmg_resources', 'version']) / 'parsed_spec.json'}

    def queue_jobs(
        self,
        multicohort: targets.MultiCohort,
        inputs: stage.StageInput,
    ) -> stage.StageOutput:
        output = self.expected_outputs(multicohort)

        batch_instance = hail_batch.get_batch()

        mane_input = batch_instance.read_input(config.config_retrieve(['acmg_resources', 'mane_summary']))

        specification_input = batch_instance.read_input(config.config_retrieve(['acmg_resources', 'specification']))

        job = batch_instance.new_bash_job(name='Parse ACMG Specification', attributes=self.get_job_attrs(multicohort))
        job.image(config.config_retrieve(['workflow', 'driver_image']))
        job.command(f"""
        python -m talos_af.scripts.process_acmg_spec \\
            --input {specification_input} \\
            --mane {mane_input} \\
            --output {job.output}
        """)

        batch_instance.write_output(job.output, output['json'])

        return self.make_outputs(multicohort, output, jobs=job)


@stage.stage(required_stages=ParseAcmgSpec)
class GenerateBedFromAcmg(stage.MultiCohortStage):
    """
    Use the ACMG specification and a GFF3 file to generate a BED file, representing this spec's region of interest.
    The NF workflow includes a fresh generation of this and a region filter, but this is required to export VCF from MT
    """

    def expected_outputs(self, dataset: targets.MultiCohort) -> dict[str, Path]:
        return {'bed': self.prefix / config.config_retrieve(['acmg_resources', 'version']) / 'parsed_spec.bed'}

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        batch_instance = hail_batch.get_batch()

        output = self.expected_outputs(multicohort)

        local_gff3 = batch_instance.read_input(config.config_retrieve(['acmg_resources', 'mane_summary']))

        parsed_spec = batch_instance.read_input(inputs.as_str(multicohort, ParseAcmgSpec))

        batch_instance = hail_batch.get_batch()

        job = batch_instance.new_bash_job(name='Parse GFF3 into BED', attributes=self.get_job_attrs(multicohort))
        job.image(config.config_retrieve(['workflow', 'driver_image']))
        job.command(f"""
        python -m talos_af.scripts.process_gff3 \
            --input {parsed_spec} \
            --gff3 {local_gff3} \
            --output {job.output}
        """)
        batch_instance.write_output(job.output, output['bed'])

        return self.make_outputs(multicohort, output, jobs=job)


@stage.stage(required_stages=GenerateBedFromAcmg)
class ExportVcfFromMt(stage.DatasetStage):
    """Find the latest AnnotateDataset output for this Dataset, export it as a VCF."""

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        # dataset_bucket / wf_name / ACMG_version / stage_name / vcf
        return {
            'vcf': dataset.prefix()
            / workflow.get_workflow().name
            / config.config_retrieve(['acmg_resources', 'version'])
            / self.name
            / 'region_filtered.vcf.bgz'
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(dataset)

        batch_instance = hail_batch.get_batch()

        input_mt = internal_utils.query_for_latest_analysis(dataset=dataset.name, stage_name='AnnotateDataset')

        bed_file = inputs.as_str(dataset, GenerateBedFromAcmg)

        job = batch_instance.new_bash_job(f'VCF from MT: {dataset.name}', attributes=self.get_job_attrs(dataset))
        job.image(config.config_retrieve(['workflow', 'driver_image']))
        job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
        job.command(f"""
            python -m talos_af.scripts.extract_vcf_from_mt \\
            --input {input_mt} \\
            --output {job.output} \\
            --bed {bed_file}
        """)

        batch_instance.write_output(job.output, str(output['vcf']).removesuffix('.vcf.bgz'))

        return self.make_outputs(dataset, output, jobs=job)


def main(dry_run: bool = False):
    if not dry_run:
        hail_batch.get_batch(attributes={'talos_af': 'true'})

    workflow.run_workflow(
        name='talos_af',
        stages=[ExportVcfFromMt],
        dry_run=dry_run,
    )


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--dry-run', action='store_true', help='Print the commands that would be run')
    args = parser.parse_args()
    main(dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
