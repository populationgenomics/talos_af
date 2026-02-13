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

    # region: download clinvar source files
    # this could be written more elegantly, but it's not worth thinking through
    submission_name = 'submission_summary.txt.gz'
    if (clinvar_file_tmp / submission_name).exists():
        job.subs = batch_instance.read_input(clinvar_file_tmp / submission_name)
    else:
        job.command(f'wget -q {CLINVAR_FTP}/{submission_name} -O {job.subs}')
        batch_instance.write_output(job.subs, clinvar_file_tmp / submission_name)

    variant_name = 'variant_summary.txt.gz'
    if (clinvar_file_tmp / variant_name).exists():
        job.vars = batch_instance.read_input(clinvar_file_tmp / variant_name)
    else:
        job.command(f'wget -q {CLINVAR_FTP}/{variant_name} -O {job.vars}')
        batch_instance.write_output(job.vars, clinvar_file_tmp / variant_name)
    # endregion

    # region: make new summary
    job.command(f"""
        python3 -m clinvarbitration.scripts.resummarise_clinvar \\
            -v {job.vars} \\
            -s {job.subs} \\
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
