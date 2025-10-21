process MakeSitesOnlyVcf {
    container params.container

    // take the VCF and index
    input:
        tuple path(vcf), path(tbi)

    // strip this VCF down to sites-only
    publishDir params.output_dir

    output:
        tuple \
            path("${params.cohort}_sites_only.vcf.bgz"), \
            path("${params.cohort}_sites_only.vcf.bgz.tbi")

    script:
    """
    bcftools view \
        --no-version \
        -W=tbi \
        -G \
        -Oz -o "${params.cohort}_sites_only.vcf.bgz" \
        ${vcf}
    """
}
