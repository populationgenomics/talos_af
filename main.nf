#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AnnotateCsq } from './nextflow/modules/talos_af/AnnotateCsq/main'
include { AnnotateWithEchtvar } from './nextflow/modules/talos_af/AnnotateWithEchtvar/main'
include { ApplyMoiFiltering } from './nextflow/modules/talos_af/ApplyMoiFiltering/main'
include { FilterVcfToBed } from './nextflow/modules/talos_af/FilterVcfToBed/main'
include { GenerateHtmlReport } from './nextflow/modules/talos_af/GenerateHtmlReport/main'

def timestamp = new java.util.Date().format('yyyy-MM')

workflow {
    main:

    if (file(workflow.outputDir).simpleName == file(params.processed_annotations).simpleName) {
    	println "Output Directory (${workflow.outputDir}) is probably not set correctly, use config or `-output-dir`"
		exit 1
    }

    println "${workflow.outputDir}"
    println "${params.processed_annotations}"

    if (file(params.acmg_bed).exists()) {
        ch_acmg_bed = Channel.fromPath(params.acmg_bed, checkIfExists: true)
		ch_acmg_json = Channel.fromPath(params.acmg_json, checkIfExists: true)
    }
    else {
    	println "ACMG artifacts (BED, JSON) don't exist, please run the Prep workflow (preparation.nf)"
		exit 1
    }

	ch_ref_genome = Channel.fromPath(params.ref_genome, checkIfExists: true)
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

	ch_alphamissense_echtvar = Channel.fromPath(params.alphamissense_echtvar)

    // does this month's clinvarbitration data exist?
    String current_clinvarbitration = "${params.processed_annotations}/clinvarbitration_${timestamp}.zip"

    if (file(current_clinvarbitration).exists()) {
        ch_clinvar_echtvar = Channel.fromPath(current_clinvarbitration)
    }
    else {
    	println "ClinvArbitration data for the month is not available, please run the Prep workflow (preparation.nf)"
		exit 1
    }

	ch_revel_echtvar = Channel.fromPath(params.revel_echtvar)

    FilterVcfToBed(
        ch_input_vcf,
        ch_input_vcf_index,
        ch_acmg_bed,
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
        ch_acmg_json,
        ch_config,
        ch_previous_results,
    )

    GenerateHtmlReport(
        ApplyMoiFiltering.out,
        ch_config,
    )

    publish:
    	vcf = AnnotateCsq.out
    	report = GenerateHtmlReport.out.report
}

output {
    vcf {
    }
    report {
    }
}
