process AnnotateCsq {
    container params.container

    input:
        path vcf
        path gff3
        path reference

    // annotate this VCF with gnomAD data
    publishDir params.output_dir, mode: 'copy'

    output:
        tuple \
            path("${params.cohort}_csq.vcf.bgz"), \
            path("${params.cohort}_csq.vcf.bgz.tbi")

    script:
    """
    bcftools index -t ${vcf}
    bcftools csq \
        --no-version \
        --force \
        -f "${reference}" \
        --local-csq \
        -g ${gff3} \
        -B 20 \
        -W=tbi \
        -Oz -o "${params.cohort}_csq.vcf.bgz" \
        ${vcf}
    """
}
