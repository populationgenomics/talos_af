process GenerateHtmlReport {
    container params.container

    input:
        path results
        path config
        path ext_ids

    output:
        path "${params.cohort}_results.html", emit: report

    script:

	def igv_arg = params.igv_dir != 'UNDEFINED' ? "--igv_dir ${params.igv_dir}" : ''
	def ext_id_arg = ext_ids.name != 'NO_FILE' ? "--ext_ids ${ext_ids}" : ''

        """
        export TALOS_AF_CONFIG=${config}

        python -m talos_af.generate_report \
            --input ${results} \
            --output ${params.cohort}_results.html $igv_arg $ext_id_arg
        """
}
