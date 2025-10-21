process ParseGff3IntoBed {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
        path acmg
        path gff3

    // generates a TSV and HT representation of the contents
    output:
        path "acmg_regions.bed"

    script:
        """
        python -m talos_af.scripts.process_gff3 \
            --input ${acmg} \
            --gff3 ${gff3} \
            --output acmg_regions.bed
        """
}
