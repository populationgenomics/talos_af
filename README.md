# Talos-AF

A [Talos](https://github.com/populationgenomics/talos)-like application for applying ACMG Actionable Finding detection to large cohorts

Removing Hail, hopefully not necessary - this is a really small region of interest, so the throughput optimisations of
Hail shouldn't be worth the effort of adding that as a dependency

1. parse the ACMG specification, producing a JSON file of each gene, and the corresponding interpretation
2. parse the GFF3 file for the genome build, and filter to genes of interest. Make a BED file of genes +/- a span
3. taking an input (joint-called) VCF, filter to the BED regions
4. annotate with gnomAD using echtvar
5. annotate with clinvar using bcftools annotate
6. annotate transcript consequences using bcftools csq
7. ...
8. PROFIT!
