process PrepareAcmgSpec {
    container params.container

    input:
        path spec
		path mane

    // read the ACMG specification, and generate a JSON summary
    output:
        path "parsed_acmg.json", emit: json
        path "parsed_acmg.bed", emit: bed

    script:
        """
        python -m talos_af.scripts.process_acmg_spec \
            --input ${spec} \
            --mane ${mane} \
            --json_out parsed_acmg.json \
            --bed_out parsed_acmg.bed
        """
}
