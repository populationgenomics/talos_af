process EncodeClinvar {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path vcf
        path encode_json

    output:
        path "filtered_clinvar.zip"

    script:
        """
        echtvar encode filtered_clinvar.zip ${encode_json} ${vcf}
        """
}
