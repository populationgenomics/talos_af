process EncodeClinvar {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path vcf

    output:
        path "filtered_clinvar.zip"

    script:
        """
        echtvar encode filtered_clinvar.zip /talos_af/echtvar/clinvar_config.json ${vcf}
        """
}
