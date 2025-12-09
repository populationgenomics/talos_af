process PrepareAcmgSpec {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
        path spec
		path mane

    // read the ACMG specification, and generate a JSON summary
    output:
        path "parsed_acmg.json"

    script:
        """
        python -m talos_af.scripts.process_acmg_spec \
            --input ${spec} \
            --mane ${mane} \
            --output parsed_acmg.json
        """
}
