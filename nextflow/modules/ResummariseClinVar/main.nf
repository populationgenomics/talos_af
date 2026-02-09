process ResummariseClinVar {
    container params.container

    input:
        path variant_summary
        path submission_summary
        val timestamp

    output:
        path "clinvarbitration_${timestamp}.vcf.bgz", emit: "vcf"
        path "clinvarbitration_${timestamp}.vcf.bgz.tbi", emit: "vcf_idx"

    script:
    """
    python3 -m clinvarbitration.scripts.resummarise_clinvar \
        -v "${variant_summary}" \
        -s "${submission_summary}" \
        -o "clinvarbitration_wasted" \
        --all_vcf "clinvarbitration_${timestamp}.vcf.bgz"

    # remove all the byproducts I don't need
    rm -r clinvarbitration_wasted*
    """
}
