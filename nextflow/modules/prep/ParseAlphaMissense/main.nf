process ParseAlphaMissense {
    container params.container

    input:
        path am_tsv
        path bed

    output:
        path "filtered_alphamissense.vcf.gz"

    script:
        """
        python -m talos_af.scripts.process_alphamissense \
            --input ${am_tsv} \
            --regions ${bed} \
            --output filtered_alphamissense.vcf.gz
        """
}
