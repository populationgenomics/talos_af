#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AnnotateCsq } from './nextflow/modules/AnnotateCsq/main'
include { AnnotateWithEchtvar } from './nextflow/modules/AnnotateWithEchtvar/main'
include { ApplyMoiFiltering } from './nextflow/modules/ApplyMoiFiltering/main'
include { EncodeAlphaMissense } from './nextflow/modules/EncodeAlphaMissense/main'
include { EncodeRevel } from './nextflow/modules/EncodeRevel/main'
include { FilterVcfToBed } from './nextflow/modules/FilterVcfToBed/main'
include { GenerateHtmlReport } from './nextflow/modules/GenerateHtmlReport/main'
include { ParseAlphaMissense } from './nextflow/modules/ParseAlphaMissense/main'
include { ParseRevel } from './nextflow/modules/ParseRevel/main'
include { PrepareAcmgSpec } from './nextflow/modules/PrepareAcmgSpec/main'

// download and prepare ClinVar data
include { DownloadClinVarFiles } from './nextflow/modules/DownloadClinVarFiles/main'
include { EncodeClinvar } from './nextflow/modules/EncodeClinvar/main'
include { ResummariseClinVar } from './nextflow/modules/ResummariseClinVar/main'

workflow {
    main:

    def timestamp = new java.util.Date().format('yyyy-MM')

    // populate various input channels
    ch_ref_genome = Channel.fromPath(params.ref_genome, checkIfExists: true)
	ch_acmg_spec = Channel.fromPath(params.acmg_spec, checkIfExists: true)
	ch_mane_input = Channel.fromPath(params.mane_input, checkIfExists: true)
	ch_gff3 = Channel.fromPath(params.gff_input, checkIfExists: true)

    ch_input_vcf = Channel.fromPath(params.input_vcf, checkIfExists: true)
    ch_input_vcf_index = Channel.fromPath("${params.input_vcf}.tbi", checkIfExists: true)
    ch_pedigree = Channel.fromPath(params.pedigree, checkIfExists: true)
    ch_config = Channel.fromPath(params.config, checkIfExists: true)

    // read the echtvar reference file as an input channel
    ch_gnomad_echtvar = Channel.fromPath(params.gnomad_echtvar, checkIfExists: true)

    // optional path to a previous set of results
    ch_previous_results = Channel.fromPath(params.previous_results, checkIfExists: true)

    PrepareAcmgSpec(
        ch_acmg_spec,
        ch_mane_input,
    )

    // generate the AlphaMissense zip if it doesn't already exist
    if (file(params.alphamissense_echtvar).exists()) {
        ch_alphamissense_echtvar = Channel.fromPath(params.alphamissense_echtvar)
    }
    else {
    	ch_alphamissense_input = Channel.fromPath(params.alphamissense_input, checkIfExists: true)
        ParseAlphaMissense(
            ch_alphamissense_input,
            PrepareAcmgSpec.out.bed,
        )
        EncodeAlphaMissense(ParseAlphaMissense.out)
        ch_alphamissense_echtvar = EncodeAlphaMissense.out
    }

    // generate the Clinvar zip if it doesn't already exist

    // does this month's clinvarbitration data exist?
    String current_clinvarbitration = "${params.processed_annotations}/clinvarbitration_${timestamp}.zip"

    if (file(params.clinvar_echtvar).exists()) {
        ch_clinvar_echtvar = Channel.fromPath(params.clinvar_echtvar)
    }
    else {
    // new workflow elements to go and create it from raw data
        String subfile = "${params.large_files}/submissions_${timestamp}.txt.gz"
        String varfile = "${params.large_files}/variants_${timestamp}.txt.gz"

        if (file(subfile).exists() && file(varfile).exists()) {
            ch_clinvar_sub = Channel.fromPath(subfile)
            ch_clinvar_var = Channel.fromPath(varfile)
        } else {
            DownloadClinVarFiles(timestamp)
            ch_clinvar_sub = DownloadClinVarFiles.out.submissions
            ch_clinvar_var = DownloadClinVarFiles.out.variants
        }

        ResummariseClinVar(
            ch_clinvar_var,
            ch_clinvar_sub,
            timestamp,
        )
        EncodeClinvar(
            ResummariseClinVar.out.vcf,
        )
        ch_clinvar_echtvar = EncodeClinvar.out
    }

    // generate the REVEL zip if it doesn't already exist
    if (file(params.revel_echtvar).exists()) {
        ch_revel_echtvar = Channel.fromPath(params.revel_echtvar)
    }
    else {
    	ch_revel_input = Channel.fromPath(params.revel_input, checkIfExists: true)
        ParseRevel(
            ch_revel_input,
            PrepareAcmgSpec.out.bed,
        )
        EncodeRevel(ParseRevel.out)
        ch_revel_echtvar = EncodeRevel.out
    }

    FilterVcfToBed(
        ch_input_vcf,
        ch_input_vcf_index,
        PrepareAcmgSpec.out.bed,
        ch_ref_genome,
    )

    AnnotateWithEchtvar(
        FilterVcfToBed.out,
        ch_gnomad_echtvar,
        ch_revel_echtvar,
        ch_alphamissense_echtvar,
        ch_clinvar_echtvar,
    )

    // annotate transcript consequences with bcftools csq
    AnnotateCsq(
        AnnotateWithEchtvar.out,
        ch_gff3,
        ch_ref_genome,
    )

    // apply per-gene rules
    ApplyMoiFiltering(
        AnnotateCsq.out,
        ch_pedigree,
        PrepareAcmgSpec.out.json,
        ch_config,
        ch_previous_results,
    )

    GenerateHtmlReport(
        ApplyMoiFiltering.out,
        ch_config,
    )
}
