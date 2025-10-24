process ParseClinvar {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path clinvar_tar
        path bed
        path header

    output:
        path "filtered_clinvar.vcf.bgz"

    script:
        """
        tar --no-same-owner -zxf ${clinvar_tar}

        python -m talos_af.scripts.process_clinvar \
            --input clinvarbitration_data/clinvar_decisions.tsv \
            --regions ${bed} \
            --header ${header} |
        bgzip -c > filtered_clinvar.vcf.bgz
        tabix filtered_clinvar.vcf.bgz
        """
}
