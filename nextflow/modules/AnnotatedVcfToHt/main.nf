process AnnotatedVcfToHt {
    container params.container

    input:
        tuple path(vcf), path(tbi)
        path alphamissense
        path revel
        path acmg_spec

    publishDir params.output_dir

    output:
        path "${params.cohort}_annotations.ht"

    script:
        """
        set -ex

        python -m talos_af.scripts.process_sitesonly_vcf \
            --input ${vcf} \
            --acmg_spec ${acmg_spec} \
            --am ${alphamissense} \
            --revel ${revel} \
            --output ${params.cohort}_annotations.ht
        """
}
