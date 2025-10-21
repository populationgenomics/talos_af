
process AnnotateGnomadAfs {
    container params.container

    input:
        tuple path(vcf), path(tbi)
        path zip

    // annotate this VCF with gnomAD data
    publishDir params.output_dir

    output:
        path "${params.cohort}_gnomad.vcf.bgz"

    script:
        """
        set -ex
        echtvar anno \
            -e ${zip} \
            ${vcf} \
            "${params.cohort}_gnomad.vcf.bgz"
        """
}
