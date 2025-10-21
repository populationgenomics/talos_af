process CombineAnnotationsAndFilterMt {
    container params.container

    input:
        tuple path(vcf), path(tbi)
        path annotations
        path clinvar_tar
        path acmg_spec


    publishDir params.output_dir, mode: 'copy'

    output:
        path "${params.cohort}.mt"

    script:
        """
        set -ex

        tar --no-same-owner -zxf ${clinvar}

        python -m talos_af.scripts.process_annotated_callset \
            --input ${vcf} \
            --annotations ${annotations} \
            --clinvar clinvarbitration_data/clinvar_decisions.ht \
            --acmg_spec ${acmg_spec} \
            --output ${params.cohort}.mt
        """
}
