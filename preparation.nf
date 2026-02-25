#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EncodeAlphaMissense } from './nextflow/modules/prep/EncodeAlphaMissense/main'
include { EncodeRevel } from './nextflow/modules/prep/EncodeRevel/main'
include { ParseAlphaMissense } from './nextflow/modules/prep/ParseAlphaMissense/main'
include { ParseRevel } from './nextflow/modules/prep/ParseRevel/main'
include { PrepareAcmgSpec } from './nextflow/modules/prep/PrepareAcmgSpec/main'

// download and prepare ClinVar data
include { DownloadClinVarFiles } from './nextflow/modules/DownloadClinVarFiles/main'
include { EncodeClinvar } from './nextflow/modules/EncodeClinvar/main'
include { ResummariseClinVar } from './nextflow/modules/ResummariseClinVar/main'

def timestamp = new java.util.Date().format('yyyy-MM')

workflow {
    main:
	ch_acmg_spec = Channel.fromPath(params.acmg_spec, checkIfExists: true)
	ch_mane_input = Channel.fromPath(params.mane_input, checkIfExists: true)

    PrepareAcmgSpec(
        ch_acmg_spec,
        ch_mane_input,
    )

	ch_alphamissense_input = Channel.fromPath(params.alphamissense_input, checkIfExists: true)
	ParseAlphaMissense(
		ch_alphamissense_input,
		PrepareAcmgSpec.out.bed,
	)
	EncodeAlphaMissense(ParseAlphaMissense.out)

    // does this month's clinvarbitration data exist? Can be downloaded using large_files/download_inputs.sh
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

	ch_revel_input = Channel.fromPath(params.revel_input, checkIfExists: true)
	ParseRevel(
		ch_revel_input,
		PrepareAcmgSpec.out.bed,
	)
	EncodeRevel(ParseRevel.out)

    // use workflow outputs, not individual copies
    publish:
    	alphamissense = EncodeAlphaMissense.out
    	clinvar = EncodeClinvar.out
    	clinvar_sub = ch_clinvar_sub
    	clinvar_var = ch_clinvar_var
        revel = EncodeRevel.out
}

output {
    alphamissense {
    }
    clinvar {
    }
    clinvar_sub {
    }
    clinvar_var {
    }
    revel {
    }
}
