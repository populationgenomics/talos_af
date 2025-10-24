process EncodeAlphaMissense {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path vcf
        path encode_json

    output:
        path "filtered_alphamissense.zip"

    script:
        """
        echtvar encode filtered_alphamissense.zip ${encode_json} ${vcf}
        """
}
