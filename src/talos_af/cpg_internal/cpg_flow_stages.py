from argparse import ArgumentParser

from cpg_flow import stage, targets, workflow
from cpg_flow.stage import StageInput, StageOutput
from cpg_flow.targets import MultiCohort
from cpg_utils import Path, config, hail_batch, to_path

from talos_af import utils as utils_af
from talos_af.cpg_internal import utils as utils_internal

ACMG_VERSION = config.config_retrieve(['acmg_resources', 'version'])


@stage.stage
class ParseAcmgSpec(stage.MultiCohortStage):
    """
    Read the TSV representing the ACMG specification, and parse it
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        return {'json': self.prefix / ACMG_VERSION / 'parsed_spec.json'}

    def queue_jobs(
        self,
        multicohort: targets.MultiCohort,
        inputs: stage.StageInput,
    ) -> stage.StageOutput:
        output = self.expected_outputs(multicohort)

        batch_instance = hail_batch.get_batch()

        mane_input = batch_instance.read_input(config.config_retrieve(['references', 'mane_summary']))

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
        return {'bed': self.prefix / ACMG_VERSION / 'parsed_spec.bed'}

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        batch_instance = hail_batch.get_batch()

        output = self.expected_outputs(multicohort)

        local_gff3 = batch_instance.read_input(config.config_retrieve(['references', 'gff3']))

        parsed_spec = batch_instance.read_input(inputs.as_str(multicohort, ParseAcmgSpec, 'json'))

        batch_instance = hail_batch.get_batch()

        job = batch_instance.new_bash_job(name='Parse GFF3 into BED', attributes=self.get_job_attrs(multicohort))
        job.image(config.config_retrieve(['workflow', 'driver_image']))
        job.command(f"""
        python -m talos_af.scripts.process_gff3 \\
            --input {parsed_spec} \\
            --gff3 {local_gff3} \\
            --output {job.output}
        """)
        batch_instance.write_output(job.output, output['bed'])

        return self.make_outputs(multicohort, output, jobs=job)


@stage.stage(required_stages=GenerateBedFromAcmg)
class GenerateRevelZip(stage.MultiCohortStage):
    def expected_outputs(self, dataset: targets.MultiCohort) -> dict[str, Path]:
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

        bed_local = batch_instance.read_input(inputs.as_str(multicohort, GenerateBedFromAcmg, 'bed'))

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


@stage.stage(required_stages=GenerateBedFromAcmg)
class GenerateAlphaMissenseZip(stage.MultiCohortStage):
    def expected_outputs(self, dataset: targets.MultiCohort) -> dict[str, Path]:
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

        bed_local = batch_instance.read_input(inputs.as_str(multicohort, GenerateBedFromAcmg, 'bed'))

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
                --output filtered_revel.vcf.gz
            """)

        # and encode that result as an Echtvar resource
        job.command(f"""
            echtvar encode {job.output['echtvar.zip']} /talos_af/echtvar/am_config.json filtered_revel.vcf.gz
        """)

        batch_instance.write_output(job.output, str(outputs['raw']).removesuffix('.raw.zip'))

        return self.make_outputs(multicohort, outputs, jobs=job)


@stage.stage(required_stages=GenerateBedFromAcmg)
class GenerateClinvarZip(stage.MultiCohortStage):
    def expected_outputs(self, dataset: targets.MultiCohort) -> dict[str, Path]:
        zenodo_link = config.config_retrieve(['references', 'clinvar_record'])
        rec_id, _dl = utils_af.get_latest_zenodo_file(zenodo_link)
        prefix = to_path(config.dataset_path(dataset='common', suffix='references/acmg_actionable')) / ACMG_VERSION
        return {
            'raw': prefix / f'{rec_id}.tar.gz',
            'zip': prefix / f'{rec_id}.echtvar.zip',
        }

    def queue_jobs(
        self,
        multicohort: MultiCohort,
        inputs: StageInput,
    ) -> StageOutput:
        outputs = self.expected_outputs(multicohort)

        batch_instance = hail_batch.get_batch()

        _id, dl = utils_af.get_latest_zenodo_file(config.config_retrieve(['references', 'clinvar_record']))

        bed_local = batch_instance.read_input(inputs.as_str(multicohort, GenerateBedFromAcmg, 'bed'))

        job = batch_instance.new_bash_job('Download and process ClinvArbitration data')
        job.image(config.config_retrieve(['workflow', 'driver_image']))

        job.declare_resource_group(
            output={
                'tar.gz': '{root}.tar.gz',
                'echtvar.zip': '{root}.echtvar.zip',
            }
        )

        # download the raw tar file
        job.command(f'curl -o {job.output["tar.gz"]} {dl}')

        # convert decisions to a VCF, and region filter
        job.command(f"""
            tar --no-same-owner -zxf {job.output['tar.gz']}

            python -m talos_af.scripts.process_clinvar \\
                --input clinvarbitration_data/clinvar_decisions.tsv \\
                --regions {bed_local} \\
                --output filtered_clinvar.vcf.gz
            """)

        # and encode that result as an Echtvar resource
        job.command(f"""
            echtvar encode {job.output['echtvar.zip']} /talos_af/echtvar/clinvar_config.json filtered_clinvar.vcf.gz
        """)

        batch_instance.write_output(job.output, str(outputs['raw']).removesuffix('.tar.gz'))

        return self.make_outputs(multicohort, outputs, jobs=job)


@stage.stage
class ExportMtFromVds(stage.DatasetStage):
    """Optional, generate a starting MT from a VDS."""

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        if vds := config.config_retrieve(['workflow', 'use_vds']):
            return {'mt': dataset.prefix() / workflow.get_workflow().name / self.name / f'{to_path(vds).name}.mt'}
        return {}

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(dataset)

        if not (vds := config.config_retrieve(['workflow', 'use_vds'])):
            return self.make_outputs(dataset, output)

        batch_instance = hail_batch.get_batch()

        job = batch_instance.new_bash_job(f'MT from VDS: {dataset.name}', attributes=self.get_job_attrs(dataset))
        job.image(config.config_retrieve(['workflow', 'driver_image']))
        job.command(f'python -m talos_af.scripts.mt_from_vds --input {vds} --output {output["mt"]}')
        return self.make_outputs(dataset, output, jobs=job)


@stage.stage(required_stages=[GenerateBedFromAcmg, ExportMtFromVds])
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

        if config.config_retrieve(['workflow', 'use_vds']):
            input_mt = inputs.as_str(dataset, ExportMtFromVds, 'mt')
        else:
            input_mt = utils_internal.query_for_latest_analysis(dataset=dataset.name, stage_name='AnnotateDataset')

        bed_file = inputs.as_str(workflow.get_multicohort(), GenerateBedFromAcmg, 'bed')

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
        return {
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

        # nextflow go brrrr
        job.command(
            f"""
            nextflow \
                -c nextflow/talos_af.config \\
                run nextflow/talos_af.nf \\
                --config nextflow/inputs/config.toml \\
                --pedigree {pedigree} \\
                --input_vcf {vcf_with_index} \\
                --acmg_spec {acmg_spec} \\
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

        # set some resource params
        job.storage('100Gi').memory('highmem').cpu(2)

        # copy the outputs back, in one smooooooth motion
        batch_instance.write_output(job.output, str(outputs['json']).removesuffix('_results.json'))
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
