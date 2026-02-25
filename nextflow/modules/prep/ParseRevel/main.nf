process ParseRevel {
    container params.container

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
