#!/usr/bin/env nextflow

/*
This workflow is a minimised annotation process for the Talos pipeline.

It takes multiple VCFs, merges them into a single VCF, and annotates the merged VCF with relevant annotations.

The specific annotations are:
- gnomAD v4.1 frequencies, applied to the joint VCF using echtvar
- Transcript consequences, using BCFtools annotate
- AlphaMissense, applied using Hail
- MANE trancript IDs and corresponding proteins, applied using Hail
*/

nextflow.enable.dsl=2

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
}
