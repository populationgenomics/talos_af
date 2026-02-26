process EncodeRevel {
    container params.container

    input:
        path vcf

    output:
        path "filtered_revel.zip"

    script:
        """
        echtvar encode filtered_revel.zip /talos_af/echtvar/revel_config.json ${vcf}
        """
}
