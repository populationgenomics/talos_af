from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


CLINVAR_FTP = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited'


def run_clinvarbitration_in_full(clinvar_file_tmp: Path, final_output: Path) -> 'BashJob':
    """
    gets the remote resources for submissions and variants

    Args:
        clinvar_file_tmp (Pathlike): temp directory for clinvar files
        final_output (Pathlike): where to write the encoded zip
    """

    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_bash_job('Run ClinvArbitration').storage('10Gi').memory('highmem').cpu('4')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command('set -eo pipefail')

    submission_name = 'submission_summary.txt.gz'
    variant_name = 'variant_summary.txt.gz'

    # region: download clinvar source files
    # for each input, if the file exists in GCP temp, use it, otherwise download a fresh one, write it to temp and GCP
    for filename in [submission_name, variant_name]:
        temp_gcp_file = clinvar_file_tmp / filename
        if not temp_gcp_file.exists():
            # download with wget, write to stdout, use tee to write that to GCP temp, and a local file. Nice.
            job.command(
                f'wget -q {CLINVAR_FTP}/{filename} -O {job[filename]}',
            )
            batch_instance.write_output(job[filename], temp_gcp_file)
        else:
            job[filename] = batch_instance.read_input(temp_gcp_file)
    # endregion

    # region: make new summary
    job.command(f"""
        python3 -m clinvarbitration.scripts.resummarise_clinvar \\
            -v {job[variant_name]} \\
            -s {job[submission_name]} \\
            -o $BATCH_TMPDIR/clinvarbitration \\
            --all_vcf $BATCH_TMPDIR/clinvarbitration.vcf.bgz
    """)
    # endregion

    # region: encode as an echtvar zip
    job.command(
        f'echtvar encode {job.zip} /talos_af/echtvar/clinvar_config.json $BATCH_TMPDIR/clinvarbitration.vcf.bgz'
    )
    batch_instance.write_output(job.zip, final_output)
    # endregion

    return job
