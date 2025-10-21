process ParseRevelIntoHt {
    container params.container

    publishDir params.processed_annotations, mode: 'copy'

    input:
        path revel
        path bed

    // generates a TSV and HT representation of the contents
    output:
        path "filtered_revel.tsv", emit: 'tsv'
        path "filtered_revel.ht", emit: 'ht'

    script:
        """
        python -m talos_af.scripts.process_revel \
            --input ${revel} \
            --bed ${bed} \
            --output_tsv filtered_revel.tsv \
            --output_ht filtered_revel.ht
        """
}
