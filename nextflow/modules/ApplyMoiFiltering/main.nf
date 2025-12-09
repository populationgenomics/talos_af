process ApplyMoiFiltering {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
        tuple path(vcf), path(vcf_tbi)
        path pedigree
        path spec
        path config
        path previous

    output:
        path "${params.cohort}_results.json"

	script:
		def history_arg = previous.name != 'NO_FILE' ? "--prior_results $previous_" : ''

    script:
        """
        export TALOS_AF_CONFIG=${config}

        python -m talos_af.apply_moi_filtering \
            --vcf ${vcf} \
            --acmg_spec ${spec} \
            --pedigree ${pedigree} \
            --output ${params.cohort}_results.json ${history_arg}
        """
}
