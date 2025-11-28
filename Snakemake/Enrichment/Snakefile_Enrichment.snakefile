import pandas as pd
import os

SAMPLES = config.get("SAMPLES", "").split(",") if config.get("SAMPLES") else []

rule all:
    input:
        expand("/cluster/projects/bhklab/procdata/Radiogenomics/genomics/gene_signatures/v2/{sample}_hallmark_signatures.csv", sample=SAMPLES)

rule filter_patients:
    input:
        dataset=lambda wildcards: config["datasets"][wildcards.sample],
        ids=lambda wildcards: config["patient_ids"][wildcards.sample]
    output:
        "/cluster/projects/bhklab/procdata/Radiogenomics/genomics/patient_ids/{sample}_filtered.csv"
    resources:
        mem_mb=2000,
        runtime=30
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/genomics/preprocessing/Data_sampler.R {input.dataset} {input.ids} {output}
        """

rule hallmark_signature_extraction:
    input:
        dataset="/cluster/projects/bhklab/procdata/Radiogenomics/genomics/patient_ids/{sample}_filtered.csv",
        gmt=lambda wildcards: config["hallmark_gmt"]
    output:
        "/cluster/projects/bhklab/procdata/Radiogenomics/genomics/gene_signatures/v2/{sample}_hallmark_signatures.csv"
    resources:
        mem_mb=8000,
        runtime=120
    shell:
        """
        module load R
        Rscript /cluster/home/t138199uhn/Scripts/genomics/enrichment/hallmark_enrichment_GMT.R {input.dataset} {input.gmt} {output}
        """

