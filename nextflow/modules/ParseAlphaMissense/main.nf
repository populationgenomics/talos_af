process ParseAlphaMissense {
    container params.container

    // parse AM data as a VCF, then encode as a condensed zip file
    publishDir params.processed_annotations, mode: 'copy'

    input:
        path am_tsv
        path bed
        path header

    output:
        path "filtered_alphamissense.vcf.gz"

    script:
        """
        python -m talos_af.scripts.process_alphamissense \
            --input ${am_tsv} \
            --regions ${bed} \
            --header ${header} \
            --output filtered_alphamissense.vcf.gz
        """
}
