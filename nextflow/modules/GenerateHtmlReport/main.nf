process GenerateHtmlReport {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
        path results

    output:
        path "${params.cohort}_results.html"

    script:
        """
        export TALOS_AF_CONFIG=${config}

        python -m talos_af.generate_report \
            --input ${results} \
            --output ${params.cohort}_results.html
        """
}
