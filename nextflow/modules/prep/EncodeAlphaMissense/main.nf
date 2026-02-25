process EncodeAlphaMissense {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path vcf

    output:
        path "filtered_alphamissense.zip"

    script:
        """
        echtvar encode filtered_alphamissense.zip /talos_af/echtvar/am_config.json ${vcf}
        """
}
