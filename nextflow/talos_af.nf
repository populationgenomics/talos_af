#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AnnotateCsq } from './modules/AnnotateCsq/main'
include { AnnotateWithEchtvar } from './modules/AnnotateWithEchtvar/main'
include { ApplyMoiFiltering } from './modules/ApplyMoiFiltering/main'
include { EncodeAlphaMissense } from './modules/EncodeAlphaMissense/main'
include { EncodeClinvar } from './modules/EncodeClinvar/main'
include { EncodeRevel } from './modules/EncodeRevel/main'
include { FilterVcfToBed } from './modules/FilterVcfToBed/main'
include { ParseAlphaMissense } from './modules/ParseAlphaMissense/main'
include { ParseClinvar } from './modules/ParseClinvar/main'
include { ParseGff3IntoBed } from './modules/ParseGff3IntoBed/main'
include { ParseRevel } from './modules/ParseRevel/main'
include { PrepareAcmgSpec } from './modules/PrepareAcmgSpec/main'

workflow {

    // populate various input channels
    ch_ref_genome = channel.fromPath(params.ref_genome, checkIfExists: true)
	ch_acmg_spec = channel.fromPath(params.acmg_spec, checkIfExists: true)
	ch_mane_input = channel.fromPath(params.mane_input, checkIfExists: true)
	ch_gff3 = channel.fromPath(params.gff_input, checkIfExists: true)

    ch_input_vcf = channel.fromPath(params.input_vcf, checkIfExists: true)
    ch_input_vcf_index = channel.fromPath("${params.input_vcf}.tbi", checkIfExists: true)
    ch_pedigree = channel.fromPath(params.pedigree, checkIfExists: true)
    ch_config = channel.fromPath(params.config, checkIfExists: true)

    // read the echtvar reference file as an input channel
    ch_gnomad_echtvar = channel.fromPath(params.gnomad_echtvar, checkIfExists: true)

    PrepareAcmgSpec(
        ch_acmg_spec,
        ch_mane_input,
    )

    ParseGff3IntoBed(
        PrepareAcmgSpec.out,
        ch_gff3,
    )

    // generate the AlphaMissense zip if it doesn't already exist
    if (file(params.alphamissense_echtvar).exists()) {
        ch_alphamissense_echtvar = channel.fromPath(params.alphamissense_echtvar)
    }
    else {
    	ch_alphamissense_input = channel.fromPath(params.alphamissense_input, checkIfExists: true)
        ParseAlphaMissense(
            ch_alphamissense_input,
            ParseGff3IntoBed.out,
        )
        EncodeAlphaMissense(ParseAlphaMissense.out)
        ch_alphamissense_echtvar = EncodeAlphaMissense.out
    }

    // generate the Clinvar zip if it doesn't already exist
    if (file(params.clinvar_echtvar).exists()) {
        ch_clinvar_echtvar = channel.fromPath(params.clinvar_echtvar)
    }
    else {
        ch_clinvar_tar = channel.fromPath(params.clinvar, checkIfExists: true)
        ParseClinvar(
            ch_clinvar_tar,
            ParseGff3IntoBed.out,
        )
        EncodeClinvar(ParseClinvar.out)
        ch_clinvar_echtvar = EncodeClinvar.out
    }

    // generate the REVEL zip if it doesn't already exist
    if (file(params.revel_echtvar).exists()) {
        ch_revel_echtvar = channel.fromPath(params.revel_echtvar)
    }
    else {
    	ch_revel_input = channel.fromPath(params.revel_input, checkIfExists: true)
        ParseRevel(
            ch_revel_input,
            ParseGff3IntoBed.out,
        )
        EncodeRevel(ParseRevel.out)
        ch_revel_echtvar = EncodeRevel.out
    }

    FilterVcfToBed(
        ch_input_vcf,
        ch_input_vcf_index,
        ParseGff3IntoBed.out,
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
        PrepareAcmgSpec.out,
        ch_config,
    )
}
