import pandas as pd
import os

SAMPLES = config.get("SAMPLES", "").split(",") if config.get("SAMPLES") else []

rule all:
    input:
        expand("/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Hallmark_scores/{sample}_hallmark_signatures.csv", sample=SAMPLES),
        expand("/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Kegg_scores/{sample}_kegg_enrichment.csv", sample=SAMPLES),
        expand("/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Reactome_scores/{sample}_reactome_enrichment.csv", sample=SAMPLES),
        expand("/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Biocarta_scores/{sample}_biocarta_enrichment.csv", sample=SAMPLES),
        expand("/cluster/projects/bhklab/procdata/Radiogenomics/gene_signatures/{sample}_compiled_enrichment.csv", sample=SAMPLES),

rule filter_patients:
    input:
        dataset=lambda wildcards: config["datasets"][wildcards.sample],
        ids=lambda wildcards: config["patient_ids"][wildcards.sample]
    output:
        "/cluster/projects/bhklab/procdata/Radiogenomics/outputs/filtered_RNAseq/{sample}_filtered.csv"
    resources:
        mem_mb=2000,
        runtime=30
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/Data_sampler.R {input.dataset} {input.ids} {output}
        """

rule hallmark_signature_extraction:
    input:
        dataset="/cluster/projects/bhklab/procdata/Radiogenomics/outputs/filtered_RNAseq/{sample}_filtered.csv",
        gmt=lambda wildcards: config["hallmark_gmt"]
    output:
        "/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Hallmark_scores/{sample}_hallmark_signatures.csv"
    resources:
        mem_mb=8000,
        runtime=120
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/hallmark_enrichment_GMT.R {input.dataset} {input.gmt} {output}
        """

rule kegg_gene_enrichment:
    input:
        dataset="/cluster/projects/bhklab/procdata/Radiogenomics/outputs/filtered_RNAseq/{sample}_filtered.csv",
        gmt=lambda wildcards: config["kegg_gmt"]
    output:
        "/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Kegg_scores/{sample}_kegg_enrichment.csv"
    resources:
        mem_mb=8000,
        runtime=120
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/kegg_enrichment_GMT.R {input.dataset} {input.gmt} {output}
        """

rule reactome_gene_enrichment:
    input:
        dataset="/cluster/projects/bhklab/procdata/Radiogenomics/outputs/filtered_RNAseq/{sample}_filtered.csv",
        gmt=lambda wildcards: config["reactome_gmt"]
    output:
        "/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Reactome_scores/{sample}_reactome_enrichment.csv"
    resources:
        mem_mb=12000,
        runtime=840   # 12 hours
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/reactome_enrichment_GMT.R {input.dataset} {input.gmt} {output}
        """

rule biocarta_gene_enrichment:
    input:
        dataset="/cluster/projects/bhklab/procdata/Radiogenomics/outputs/filtered_RNAseq/{sample}_filtered.csv",
        gmt=lambda wildcards: config["biocarta_gmt"]
    output:
        "/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Biocarta_scores/{sample}_biocarta_enrichment.csv"
    resources:
        mem_mb=10000,
        runtime=120
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/biocarta_enrichment_GMT.R {input.dataset} {input.gmt} {output}
        """

rule compile_enrichment_scores:
    input:
        hallmark="/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Hallmark_scores/{sample}_hallmark_signatures.csv",
        kegg="/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Kegg_scores/{sample}_kegg_enrichment.csv",
        reactome="/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Reactome_scores/{sample}_reactome_enrichment.csv",
        biocarta="/cluster/projects/bhklab/procdata/Radiogenomics/outputs/Biocarta_scores/{sample}_biocarta_enrichment.csv",
        ids=lambda wildcards: config["patient_ids"][wildcards.sample]
    output:
        "/cluster/projects/bhklab/procdata/Radiogenomics/gene_signatures/{sample}_compiled_enrichment.csv"
    resources:
        mem_mb=4000,
        runtime=30
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/Enrichment_compiler.R {output} {input.ids} {input.hallmark} {input.kegg} {input.reactome} {input.biocarta}
        """
