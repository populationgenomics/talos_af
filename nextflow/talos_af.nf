#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AnnotateCsq } from './modules/AnnotateCsq/main'
include { AnnotateGnomadAfs } from './modules/AnnotateGnomadAfs/main'
include { FilterVcfToBed } from './modules/FilterVcfToBed/main'
include { MakeSitesOnlyVcf } from './modules/MakeSitesOnlyVcf/main'
include { ParseAlphaMissenseIntoHt } from './modules/ParseAlphaMissenseIntoHt/main'
include { ParseGff3IntoBed } from './modules/ParseGff3IntoBed/main'
include { ParseRevelIntoHt } from './modules/ParseRevelIntoHt/main'
include { PrepareAcmgSpec } from './modules/PrepareAcmgSpec/main'

workflow {

    // populate ref genome input channel
    ch_ref_genome = channel.fromPath(
    	params.ref_genome,
    	checkIfExists: true,
	)

	ch_acmg_spec = channel.fromPath(
	    params.acmg_spec,
	    checkIfExists: true,
	)

	ch_mane_input = channel.fromPath(
	    params.mane_input,
	    checkIfExists: true,
	)

	ch_gff3 = channel.fromPath(
	    params.gff_input,
	    checkIfExists: true,
	)

    PrepareAcmgSpec(
        ch_acmg_spec,
        ch_mane_input,
    )

    ParseGff3IntoBed(
        PrepareAcmgSpec.out,
        ch_gff3,
    )

    // generate the AlphaMissense HT if it doesn't already exist
    if (file(params.alphamissense_ht).exists()) {
        ch_alphamissense_table = channel.fromPath(params.alphamissense_ht)
    }
    else {
    	ch_alphamissense_input = channel.fromPath(
    		params.alphamissense_input,
    		checkIfExists: true,
		)
        ParseAlphaMissenseIntoHt(
            ch_alphamissense_input,
            ParseGff3IntoBed.out,
        )
        ch_alphamissense_table = ParseAlphaMissenseIntoHt.out.ht
    }

    // generate the REVEL HT if it doesn't already exist
    if (file(params.revel_ht).exists()) {
        ch_revel_table = channel.fromPath(params.revel_ht)
    }
    else {
    	ch_revel_input = channel.fromPath(
    		params.revel_input,
    		checkIfExists: true,
		)
		ParseRevelIntoHt(
            ch_revel_input,
            ParseGff3IntoBed.out,
        )
        ch_revel_table = ParseRevelIntoHt.out.ht
    }

    // now the data-centric steps
    ch_input_vcf = channel.fromPath(
        params.input_vcf,
        checkIfExists: true,
    )
    ch_input_vcf_index = channel.fromPath(
        "${params.input_vcf}.tbi",
        checkIfExists: true,
    )

    // read the echtvar reference file as an input channel
    ch_gnomad_echtvar = channel.fromPath(
		params.echtvar,
		checkIfExists: true
    )
    FilterVcfToBed(
        ch_input_vcf,
        ch_input_vcf_index,
        ParseGff3IntoBed.out,
        ch_ref_genome,
    )

    // create a sites-only version of this VCF, just to pass less data around when annotating
    MakeSitesOnlyVcf(
        FilterVcfToBed.out,
    )

    AnnotateGnomadAfs(
        MakeSitesOnlyVcf.out,
        ch_gnomad_echtvar,
    )

    // annotate transcript consequences with bcftools csq
    AnnotateCsq(
        AnnotateGnomadAfs.out,
        ch_gff3,
        ch_ref_genome,
    )
}
