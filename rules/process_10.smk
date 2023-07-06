import pandas as pd
from os.path import join
from snakemake.remote import AUTO
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

FTP = FTPRemoteProvider()
import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

BASE_FP = config["BASE_FP"]
DATASET = config["DATASET"]
SAMPLES = config["samples"].keys()


rule all:
    input:
        expand(
            join(BASE_FP, DATASET, "{sample}", "anndata", "adata_norm_processed.h5ad"),
            sample=SAMPLES,
        ),
        join(BASE_FP, DATASET, "merged","anndata", "adata_merged.h5ad"),
        expand(join(BASE_FP, DATASET, "merged", "plots", "umap_{gene}.jpg"),
                            gene=["Aplnr"])


include: "./common.smk"


rule process_scanpy:
    input:
        mtx_file=get_mtx_file,
        barcodes_file=get_barcodes_file,
        features_file=get_features_file,
    output:
        adata_h5=join(BASE_FP, DATASET, "{sample}", "anndata", "adata.h5ad"),
        adata_genes=join(BASE_FP, DATASET, "{sample}", "anndata", "csvs", "var.csv"),
        adata_cells=join(BASE_FP, DATASET, "{sample}", "anndata", "csvs", "obs.csv"),
        adata_counts=join(BASE_FP, DATASET, "{sample}", "anndata", "csvs", "X.csv"),
        soupx_groups=join(BASE_FP, DATASET, "{sample}", "qc", "soupx_groups.csv"),
        p1_displot=report(
            join(BASE_FP, DATASET, "{sample}", "qc", "plots", "count_dist.jpg"),
            category="QC",
            subcategory="{sample}",
        ),
        p2_mt_counts=report(
            join(BASE_FP, DATASET, "{sample}", "qc", "plots", "mt_counts.jpg"),
            category="QC",
            subcategory="{sample}",
        ),
        p3_scatter=report(
            join(BASE_FP, DATASET, "{sample}", "qc", "plots", "count_scatter.jpg"),
            category="QC",
            subcategory="{sample}",
        ),
        p3_scatter_after=report(
            join(
                BASE_FP,
                DATASET,
                "{sample}",
                "qc",
                "plots",
                "after_qc_count_scatter.jpg",
            ),
            category="QC",
            subcategory="{sample}",
        ),
    params:
        adata_csvs=directory(join(BASE_FP, DATASET, "{sample}", "anndata", "csvs")),
    conda:
        "../envs/scanpy.yml"
    threads: 1
    resources:
        mem_mb=8192 * 2,
        time_min=30,
    script:
        "../scripts/process_anndata.py"


rule run_soupx_and_scdblfind:
    input:
        mtx_file=get_mtx_file,
        features_file=get_features_file,
        adata_genes=join(BASE_FP, DATASET, "{sample}", "anndata", "csvs", "var.csv"),
        adata_cells=join(BASE_FP, DATASET, "{sample}", "anndata", "csvs", "obs.csv"),
        adata_counts=join(BASE_FP, DATASET, "{sample}", "anndata", "csvs", "X.csv"),
        soupx_groups=join(BASE_FP, DATASET, "{sample}", "qc", "soupx_groups.csv"),
    output:
        corrected_counts=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "corrected_counts.csv"
        ),
        doublet_data=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "doublet_data.csv"
        ),
    conda:
        "../envs/R.yml"
    threads: 1
    resources:
        mem_mb=8192 * 2,
        time_min=30,
    script:
        "../scripts/doublet_soup.R"


rule normalize_scanpy:
    input:
        adata_h5=join(BASE_FP, DATASET, "{sample}", "anndata", "adata.h5ad"),
        corrected_counts=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "corrected_counts.csv"
        ),
        doublet_data=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "doublet_data.csv"
        ),
    output:
        adata_out_path=temp(
            join(BASE_FP, DATASET, "{sample}", "anndata", "adata_norm_temp.h5ad")
        ),
        scran_data_mat=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "norm_data_mat.csv"
        ),
        cluster_groups=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "cluster_groups.csv"
        ),
    conda:
        "../envs/scanpy.yml"
    threads: 1
    resources:
        mem_mb=8192 * 2,
        time_min=30,
    script:
        "../scripts/normalize.py"


rule compute_size_factors:
    input:
        scran_data_mat=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "norm_data_mat.csv"
        ),
        cluster_groups=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "cluster_groups.csv"
        ),
    output:
        size_factors=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "size_factors.csv"
        ),
    conda:
        "../envs/R.yml"
    threads: 4
    resources:
        mem_mb=8192 * 2,
        time_min=30,
    script:
        "../scripts/compute_size_factors.R"


rule scran_normalize:
    input:
        adata_out_path=join(
            BASE_FP, DATASET, "{sample}", "anndata", "adata_norm_temp.h5ad"
        ),
        size_factors=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "size_factors.csv"
        ),
    output:
        adata_out_path=join(BASE_FP, DATASET, "{sample}", "anndata", "adata_norm.h5ad"),
        corrected_counts_filtered=join(
            BASE_FP,
            DATASET,
            "{sample}",
            "anndata",
            "csvs",
            "corrected_counts_filtered.csv",
        ),
        var_filtered=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "var_filtered.csv"
        ),
    conda:
        "../envs/scanpy.yml"
    threads: 1
    resources:
        mem_mb=8192,
        time_min=9,
    script:
        "../scripts/normalize_scran.py"


rule compute_biomial_deviance:
    input:
        corrected_counts=join(
            BASE_FP,
            DATASET,
            "{sample}",
            "anndata",
            "csvs",
            "corrected_counts_filtered.csv",
        ),
        var_filtered=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "var_filtered.csv"
        ),
    output:
        deviance_data=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "deviance_table.csv"
        ),
    conda:
        "../envs/R.yml"
    threads: 1
    resources:
        mem_mb=8192,
        time_min=9,
    script:
        "../scripts/biomial_deviance.R"


rule select_variable_features:
    input:
        deviance_data=join(
            BASE_FP, DATASET, "{sample}", "anndata", "csvs", "deviance_table.csv"
        ),
        adata_h5=join(BASE_FP, DATASET, "{sample}", "anndata", "adata_norm.h5ad"),
    output:
        adata_out_path=join(
            BASE_FP, DATASET, "{sample}", "anndata", "adata_norm_processed.h5ad"
        ),
        pca=report(
            join(BASE_FP, DATASET, "{sample}", "qc", "plots", "pca.jpg"),
            category="QC",
            subcategory="{sample}",
        ),
        umap=report(
            join(BASE_FP, DATASET, "{sample}", "qc", "plots", "umap.jpg"),
            category="QC",
            subcategory="{sample}",
        ),
        umap_cluster=report(
            join(BASE_FP, DATASET, "{sample}", "qc", "plots", "umap_cluster.jpg"),
            category="QC",
            subcategory="{sample}",
        ),
    conda:
        "../envs/scanpy.yml"
    threads: 1
    resources:
        mem_mb=8192,
        time_min=9,
    script:
        "../scripts/select_features.py"


rule merge_samples:
    input:
        adata_files=expand(join(
            BASE_FP, DATASET, "{sample}", "anndata", "adata_norm_processed.h5ad"
        ),sample = SAMPLES)
    output:
        big_adata=join(BASE_FP, DATASET, "merged","anndata", "adata_merged.h5ad")
    conda: "../envs/scanpy.yml"
    threads: 1
    resources:
        mem_mb = 32768,
        time_min = 9
    script: "../scripts/merge_adata.py"

rule compare_gene:
    input:
        adata_files=expand(join(
            BASE_FP, DATASET, "{sample}", "anndata", "adata_norm_processed.h5ad"
        ),sample = SAMPLES)
    output:
        umap_plot=report(expand(join(BASE_FP, DATASET, "merged", "plots", "umap_{gene}.jpg"),
                            gene=config["umap_interest"]), category = "results")
    params:
        genes=config["umap_interest"],
        samples = list(SAMPLES)
    conda: "../envs/scanpy.yml"
    threads: 1
    resources:
        mem_mb = 16384,
        time_min = 9
    script:
        "../scripts/plot_umap.py"