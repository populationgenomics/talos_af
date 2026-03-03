# Talos-AF

A [Talos](https://github.com/populationgenomics/talos)-like application for applying ACMG Actionable Finding detection to large cohorts

This workflow contains a full enclosed annotation pipeline - VCF in, annotation in the middle, results out. VCFs used
here do not need to have annotation by other sources.

## What does it do?

1. parse an ACMG specification, producing a mapping of each gene to its corresponding interpretation advice
2. using genes named in that specification, build a region of interest BED file (genes +/- a 500bp)
3. take an input VCF (single or multisample), filter to the BED regions
4. apply the same region of interest to REVEL, ClinVar, and AlphaMissense raw data
5. annotate the input VCF with gnomAD, REVEL, ClinVar, and AlphaMissense in one pass, using [Echtvar](https://github.com/brentp/echtvar)
6. annotate transcript consequences using bcftools csq
7. iterate through all lines of the resulting annotated VCF. Throw away variants which aren't damaging (ClinVar P/LP, or a high impact transcript consequence), and assess the remaining variants against the ACMG-suggested MOI
8. write a result file of all damaging variants which match the ACMG expected AD/AR/XL mode of inheritance, and the samples where they were seen

## So you wanna run it?

1. Build a docker image
    ```bash
     docker build -t talos_af:0.0.11 .
    ```
   
2. Download large input resources (used for annotation)
    ```bash
    cd large_files
    bash download_inputs.sh
    cd ..
    ```
   
3. [install NextFlow](https://www.nextflow.io/docs/latest/install.html)

4. Run the preparation workflow
    ```bash
    nextflow \
    -c nextflow.config \
        run preparation.nf \
        --acmg_spec nextflow/inputs/acmg_secondary.tsv
    ```

5. Run TalosAF!
    ```bash
    nextflow \
    -c nextflow.config \
        run main.nf \
        --acmg_spec nextflow/inputs/acmg_secondary.tsv \
        --input_vcf nextflow/inputs/test_1.vcf.bgz \
        --pedigree nextflow/inputs/pedigree.ped \
        --output-dir <path to write outputs>
    ```

For demonstration purposes, a [workflow config file](nextflow.config) has been populated with all required parameters for a test run.

By default, this will use the test data (a trio, 12 variants, corresponding pedigree) and run the full annotation and
interpretation workflow. Outputs will be written to folders inside `./nextflow` in this repository.

The first time this runs will be slow-ish, as the downloaded reference data is reformatted to be used as fast annotation
inputs. Outputs will be generated inside the `nextflow` folder. Subsequent reruns will make use of this preprocessing.

We advise that for real analyses you change the following config parameters - these can be edited in the config file, or
altered through NextFlow CLI `--[argument] [value]` syntax:
- `cohort`: a name for the collective group of samples being processed, used to label files and create output paths
- `input_vcf`: path to an input VCF containing one or more samples
- `pedigree`: path to a pedigree for the corresponding samples
- `output_dir`: path to a folder where results will be written
- `processed_annotations`: a folder to contain all the cohort-agnostic files (reformatted annotation files)

We also advise that the `output_dir` for analysis is outside this repository. The test workflow default location is
within this directory, but for real-world work mixing results with the codebase can lead to complicated git resolution
when updating the tool.
