process EncodeRevel {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path vcf
        path encode_json

    output:
        path "filtered_revel.zip"

    script:
        """
        echtvar encode filtered_revel.zip /talos_af/echtvar/revel_config.json ${vcf}
        """
}
