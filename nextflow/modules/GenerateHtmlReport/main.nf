process GenerateHtmlReport {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
        path results
        path config

    output:
        path "${params.cohort}_results.html"

    script:

	def igv_arg = params.igv_dir != 'UNDEFINED' ? "--igv_dir ${params.igv_dir}" : ''

        """
        export TALOS_AF_CONFIG=${config}

        python -m talos_af.generate_report \
            --input ${results} \
            --output ${params.cohort}_results.html $igv_arg
        """
}
