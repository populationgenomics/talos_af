import datetime
import functools
from argparse import ArgumentParser
from os.path import join

from cpg_flow import stage, targets, workflow
from cpg_flow.stage import StageInput, StageOutput
from cpg_flow.targets import MultiCohort
from cpg_utils import Path, config, hail_batch, to_path

from talos_af.cpg_internal import run_clinvarbitration
from talos_af.cpg_internal import utils as utils_internal

ACMG_VERSION = config.config_retrieve(['acmg_resources', 'version'])

THIS_MONTH = datetime.datetime.now().strftime('%y-%m')  # noqa: DTZ005


@functools.cache
def get_clinvarbitration_folder(temp_folder: bool = False) -> Path:
    """
    get the folder to use for this run
    sits in a bucket accessible to the operating dataset, in a folder named "clinvarbitration/YY-MM"
    """
    return to_path(
        join(
            config.config_retrieve(['storage', 'default', 'tmp' if temp_folder else 'default']),
            'clinvarbitration',
            THIS_MONTH,
        ),
    )


@stage.stage
class ParseAcmgSpec(stage.MultiCohortStage):
    """
    Read the TSV representing the ACMG specification, and parse it
    """

    def expected_outputs(self, _multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {
            'json': self.prefix / ACMG_VERSION / 'parsed_spec.json',
            'bed': self.prefix / ACMG_VERSION / 'parsed_spec.bed',
        }

    def queue_jobs(
        self,
        multicohort: targets.MultiCohort,
        _inputs: stage.StageInput,
    ) -> stage.StageOutput:
        outputs = self.expected_outputs(multicohort)

        batch_instance = hail_batch.get_batch()

        mane_input = batch_instance.read_input(config.config_retrieve(['references', 'mane_summary']))

        specification_input = batch_instance.read_input(config.config_retrieve(['acmg_resources', 'specification']))

        job = batch_instance.new_bash_job(name='Parse ACMG Specification', attributes=self.get_job_attrs(multicohort))
        job.image(config.config_retrieve(['workflow', 'driver_image']))
        job.command(f"""
        python -m talos_af.scripts.process_acmg_spec \\
            --input {specification_input} \\
            --mane {mane_input} \\
            --json_out {job.json_out} \\
            --bed_out {job.bed_out}
        """)

        batch_instance.write_output(job.json_out, outputs['json'])
        batch_instance.write_output(job.bed_out, outputs['bed'])

        return self.make_outputs(multicohort, outputs, jobs=job)


@stage.stage(required_stages=ParseAcmgSpec)
class GenerateRevelZip(stage.MultiCohortStage):
    def expected_outputs(self, _dataset: targets.MultiCohort) -> dict[str, Path]:
        prefix = to_path(config.dataset_path(dataset='common', suffix='references/acmg_actionable')) / ACMG_VERSION
        return {
            'raw': prefix / 'revel.raw.zip',
            'zip': prefix / 'revel.echtvar.zip',
        }

    def queue_jobs(
        self,
        multicohort: MultiCohort,
        inputs: StageInput,
    ) -> StageOutput:
        outputs = self.expected_outputs(multicohort)
        batch_instance = hail_batch.get_batch()

        bed_local = batch_instance.read_input(inputs.as_str(multicohort, ParseAcmgSpec, 'bed'))

        job = batch_instance.new_bash_job('Download and process Revel data')
        job.image(config.config_retrieve(['workflow', 'driver_image']))

        job.declare_resource_group(
            output={
                'raw.zip': '{root}.raw.zip',
                'echtvar.zip': '{root}.echtvar.zip',
            }
        )

        # download the raw zip file
        revel_link = config.config_retrieve(['references', 'revel_zenodo'])
        job.command(f'curl -o {job.output["raw.zip"]} {revel_link}')

        # convert decisions to a VCF, and region filter
        job.command(f"""
            python -m talos_af.scripts.process_revel \\
                --input {job.output['raw.zip']} \\
                --regions {bed_local} \\
                --output filtered_revel.vcf.gz
            """)

        # and encode that result as an Echtvar resource
        job.command(f"""
            echtvar encode {job.output['echtvar.zip']} /talos_af/echtvar/revel_config.json filtered_revel.vcf.gz
        """)

        batch_instance.write_output(job.output, str(outputs['raw']).removesuffix('.raw.zip'))

        return self.make_outputs(multicohort, outputs, jobs=job)


@stage.stage(required_stages=ParseAcmgSpec)
class GenerateAlphaMissenseZip(stage.MultiCohortStage):
    def expected_outputs(self, _dataset: targets.MultiCohort) -> dict[str, Path]:
        prefix = to_path(config.dataset_path(dataset='common', suffix='references/acmg_actionable')) / ACMG_VERSION
        return {
            'raw': prefix / 'am.raw.zip',
            'zip': prefix / 'am.echtvar.zip',
        }

    def queue_jobs(
        self,
        multicohort: MultiCohort,
        inputs: StageInput,
    ) -> StageOutput:
        outputs = self.expected_outputs(multicohort)
        batch_instance = hail_batch.get_batch()

        bed_local = batch_instance.read_input(inputs.as_str(multicohort, ParseAcmgSpec, 'bed'))

        job = batch_instance.new_bash_job('Download and process AlphaMissense data')
        job.image(config.config_retrieve(['workflow', 'driver_image']))

        job.declare_resource_group(
            output={
                'raw.zip': '{root}.raw.zip',
                'echtvar.zip': '{root}.echtvar.zip',
            }
        )

        # download the raw zip file
        revel_link = config.config_retrieve(['references', 'alphamissense_zenodo'])
        job.command(f'curl -o {job.output["raw.zip"]} {revel_link}')

        # convert decisions to a VCF, and region filter
        job.command(f"""
            python -m talos_af.scripts.process_alphamissense \\
                --input {job.output['raw.zip']} \\
                --regions {bed_local} \\
                --output filtered_am.vcf.gz
            """)

        # and encode that result as an Echtvar resource
        job.command(f"""
            echtvar encode {job.output['echtvar.zip']} /talos_af/echtvar/am_config.json filtered_am.vcf.gz
        """)

        batch_instance.write_output(job.output, str(outputs['raw']).removesuffix('.raw.zip'))

        return self.make_outputs(multicohort, outputs, jobs=job)


@stage.stage()
class GenerateClinvarZip(stage.MultiCohortStage):
    def expected_outputs(self, _dataset: targets.MultiCohort) -> dict[str, Path]:
        return {
            'zip': get_clinvarbitration_folder() / 'echtvar.zip',
        }

    def queue_jobs(self, multicohort: MultiCohort, _inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(multicohort)
        tmp_files_folder = get_clinvarbitration_folder(temp_folder=True)
        job = run_clinvarbitration.run_clinvarbitration_in_full(
            clinvar_file_tmp=tmp_files_folder,
            final_output=outputs['zip'],
        )
        return self.make_outputs(multicohort, outputs, jobs=job)


@stage.stage
class ExportMtFromVds(stage.DatasetStage):
    """Optional, generate a starting MT from a VDS."""

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        if vds := config.config_retrieve(['workflow', 'use_vds', dataset.name]):
            return {'mt': dataset.prefix() / workflow.get_workflow().name / self.name / f'{to_path(vds).name}.mt'}
        return {}

    def queue_jobs(self, dataset: targets.Dataset, _inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(dataset)

        if not (vds := config.config_retrieve(['workflow', 'use_vds', dataset.name])):
            return self.make_outputs(dataset, output)

        batch_instance = hail_batch.get_batch()

        job = batch_instance.new_bash_job(f'MT from VDS: {dataset.name}', attributes=self.get_job_attrs(dataset))
        job.image(config.config_retrieve(['workflow', 'driver_image']))
        job.command(f'python -m talos_af.scripts.mt_from_vds --input {vds} --output {output["mt"]}')
        return self.make_outputs(dataset, output, jobs=job)


@stage.stage(required_stages=[ParseAcmgSpec, ExportMtFromVds])
class ExportVcfFromMt(stage.DatasetStage):
    """Find the latest AnnotateDataset output for this Dataset, export it as a VCF."""

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        return {
            'vcf': dataset.prefix()
            / workflow.get_workflow().name
            / ACMG_VERSION
            / self.name
            / 'region_filtered.vcf.bgz'
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(dataset)

        batch_instance = hail_batch.get_batch()

        if config.config_retrieve(['workflow', 'use_vds', dataset.name]):
            input_mt = inputs.as_str(dataset, ExportMtFromVds, 'mt')
        else:
            input_mt = utils_internal.query_for_latest_analysis(dataset=dataset.name, stage_name='AnnotateDataset')

        bed_file = inputs.as_str(workflow.get_multicohort(), ParseAcmgSpec, 'bed')

        job = batch_instance.new_bash_job(f'VCF from MT: {dataset.name}', attributes=self.get_job_attrs(dataset))
        job.image(config.config_retrieve(['workflow', 'driver_image']))

        job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
        job.command(f"""
            python -m talos_af.scripts.extract_vcf_from_mt \\
            --input {input_mt} \\
            --output {output['vcf']} \\
            --bed {bed_file}
        """)

        return self.make_outputs(dataset, output, jobs=job)


@stage.stage(required_stages=[ExportVcfFromMt, GenerateRevelZip, GenerateAlphaMissenseZip, GenerateClinvarZip])
class RunTalosAfNextFlow(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        output_folder = dataset.prefix() / workflow.get_workflow().name / ACMG_VERSION / self.name
        web_folder = dataset.web_prefix() / workflow.get_workflow().name / ACMG_VERSION / self.name
        return {
            'html': web_folder / f'{dataset.name}_results.html',
            'json': output_folder / f'{dataset.name}_results.json',
            'vcf': output_folder / f'{dataset.name}_filtered.vcf.bgz',
            'tbi': output_folder / f'{dataset.name}_filtered.vcf.bgz.tbi',
            'pedigree': output_folder / f'{dataset.name}.pedigree',
        }

        return {}

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        batch_instance = hail_batch.get_batch()

        job = batch_instance.new_bash_job(f'Run NF workflow in full for {dataset.name}')
        job.image(config.config_retrieve(['workflow', 'driver_image']))

        # set some output expectations
        job.declare_resource_group(
            output={
                f'{dataset.name}_results.json': f'{{root}}/{dataset.name}_results.json',
                f'{dataset.name}_filtered.vcf.bgz': f'{{root}}/{dataset.name}_filtered.vcf.bgz',
                f'{dataset.name}_filtered.vcf.bgz.tbi': f'{{root}}/{dataset.name}_filtered.vcf.bgz.tbi',
            },
        )

        # read in various input files
        ref_fa = batch_instance.read_input(config.config_retrieve(['workflow', 'ref_fa']))
        gff3_localised = batch_instance.read_input(config.config_retrieve(['references', 'gff3']))
        acmg_spec = batch_instance.read_input(config.config_retrieve(['acmg_resources', 'specification']))

        multicohort = workflow.get_multicohort()
        gnomad_zip = batch_instance.read_input(config.config_retrieve(['references', 'echtvar_gnomad']))
        am_zip = batch_instance.read_input(inputs.as_str(multicohort, GenerateAlphaMissenseZip, 'zip'))
        clinvar_zip = batch_instance.read_input(inputs.as_str(multicohort, GenerateClinvarZip, 'zip'))
        revel_zip = batch_instance.read_input(inputs.as_str(multicohort, GenerateRevelZip, 'zip'))

        mane = batch_instance.read_input(config.config_retrieve(['references', 'mane_summary']))

        vcf_path = inputs.as_str(dataset, ExportVcfFromMt, 'vcf')
        vcf_with_index = batch_instance.read_input_group(
            vcf=vcf_path,
            index=f'{vcf_path}.tbi',
        )['vcf']

        pedigree = batch_instance.read_input(dataset.write_ped_file(out_path=outputs['pedigree']))

        # create the expected path to CRAM files for this dataset
        cram_dir = dataset.prefix() / 'cram'

        # nextflow go brrrr
        job.command(
            f"""
            nextflow \
                -c nextflow.config \\
                run main.nf \\
                --config nextflow/inputs/config.toml \\
                --pedigree {pedigree} \\
                --input_vcf {vcf_with_index} \\
                --acmg_spec {acmg_spec} \\
                --igv_dir {cram_dir} \\
                --mane_input {mane} \\
                --gnomad_echtvar {gnomad_zip} \\
                --revel_echtvar {revel_zip} \\
                --clinvar_echtvar {clinvar_zip} \\
                --alphamissense_echtvar {am_zip} \\
                --cohort {dataset.name} \\
                --ref_genome {ref_fa} \\
                --output_dir {job.output} \\
                --gff_input {gff3_localised} \\
                -without-docker
            """,
        )
        job.command(f'mv {job.output}/{dataset.name}_results.html {job.html}')

        # set some resource params
        job.storage('100Gi').memory('highmem').cpu(2)

        # copy the outputs back, in one smooooooth motion
        batch_instance.write_output(job.output, str(outputs['json']).removesuffix('_results.json'))
        batch_instance.write_output(job.html, outputs['html'])
        return self.make_outputs(dataset, outputs, jobs=job)


def main(dry_run: bool = False):
    if not dry_run:
        hail_batch.get_batch(attributes={'talos_af': 'true'})

    workflow.run_workflow(
        name='talos_af',
        stages=[RunTalosAfNextFlow],
        dry_run=dry_run,
    )


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--dry-run', action='store_true', help='Print the commands that would be run')
    args = parser.parse_args()
    main(dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
