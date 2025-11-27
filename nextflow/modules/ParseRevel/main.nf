process ParseRevel {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path revel
        path bed

    output:
        path "filtered_revel.vcf.gz"

    script:
        """
        python -m talos_af.scripts.process_revel \
            --input ${revel} \
            --regions ${bed} \
            --output filtered_revel.vcf.gz
        """
}
