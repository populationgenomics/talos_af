process ParseAlphaMissenseIntoHt {
    container params.container

    // parse AM data as a Hail Table
    publishDir params.processed_annotations, mode: 'copy'

    input:
        path am_tsv
        path bed

    output:
        path "filtered_alphamissense.tsv", emit: 'tsv'
        path "filtered_alphamissense.ht", emit: 'ht'

    script:
        """
        python -m talos_af.scripts.process_alphamissense \
            --input ${am_tsv} \
            --bed ${bed} \
            --output_tsv filtered_alphamissense.tsv \
            --output_ht filtered_alphamissense.ht
        """
}
