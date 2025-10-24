process ParseRevel {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path revel
        path bed
        path header

    output:
        path "filtered_revel.vcf.gz"

    script:
        """
        python -m talos_af.scripts.process_revel \
            --input ${revel} \
            --regions ${bed} \
            --header ${header} \
            --output filtered_revel.vcf.gz
        """
}
