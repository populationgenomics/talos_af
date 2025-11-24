process AnnotateWithEchtvar {
    container params.container

    input:
        tuple path(vcf), path(tbi)
        path gnomad_zip
        path revel_zip
        path am_zip
        path clinvar_zip

    // annotate VCF with gnomAD, REVEL, and AlphaMissense data
    publishDir params.output_dir

    output:
        path "${params.cohort}_echtvar.vcf.bgz"

    script:
        """
        set -ex
        echtvar anno \
            -e ${gnomad_zip} \
            -e ${revel_zip} \
            -e ${am_zip} \
            -e ${clinvar_zip} \
            -i "gnomad_AF_joint < 0.05" \
            ${vcf} \
            "${params.cohort}_echtvar.vcf.bgz"
        """
}
