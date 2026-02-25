process AnnotateWithEchtvar {
    container params.container

    input:
        tuple path(vcf), path(tbi)
        path gnomad_zip
        path revel_zip
        path am_zip
        path clinvar_zip

    output:
        path "${params.cohort}_echtvar.vcf.bgz"

    // note - the gnomAD AF filter is set very high (7.5% in gnomAD 4.1) to account for a Hom-HFE common variant
    script:
        """
        set -ex
        echtvar anno \
            -e ${gnomad_zip} \
            -e ${revel_zip} \
            -e ${am_zip} \
            -e ${clinvar_zip} \
            -i "gnomad_AF_joint < 0.75" \
            ${vcf} \
            "${params.cohort}_echtvar.vcf.bgz"
        """
}
