process EncodeClinvar {
    container params.container

    input:
        path vcf
        val timestamp

    output:
        path "clinvarbitration_${timestamp}.zip"

    script:
        """
        echtvar encode "clinvarbitration_${timestamp}.zip" /talos_af/echtvar/clinvar_config.json ${vcf}
        """
}
