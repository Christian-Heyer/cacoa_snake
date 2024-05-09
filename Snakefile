import pandas as pd
from os.path import join
from snakemake.remote import AUTO
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()
import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider


#configfile: "configs/tabularmuris.yaml"

BASE_FP = config["BASE_FP"]
DATASET = config["DATASET"] 
HTTP = HTTPRemoteProvider()
rule all:
    input:
        expand(join(BASE_FP, "{dataset}", "report", "{dataset}" + "_cacoa.r.html"), 
            dataset = DATASET)

if config["download_link"]["f_type"] == "h5ad":
    rule get_h5ad_file:
        output:
            h5ad = temp(join(BASE_FP, "{dataset}", "{dataset}_raw." + config["download_link"]["f_type"]))
        params:
            link = config["download_link"]["link"]
        shell:
            """
            wget -O {output.h5ad} {params.link}
            """
    rule anndata_to_loom:
        input:
            h5ad = join(BASE_FP, "{dataset}", "{dataset}_raw." + config["download_link"]["f_type"])
        output:
            #temp(dir(join(BASE_FP, "{dataset}", "{dataset}_raw"))), 
            counts = join(BASE_FP, "{dataset}", "{dataset}_raw", "X.csv"),
            cell_metadata = join(BASE_FP, "{dataset}", "{dataset}_raw", "obs.csv"),
            gene_metadata = join(BASE_FP, "{dataset}", "{dataset}_raw", "var.csv"),
            h5ad = join(BASE_FP, "{dataset}", "{dataset}_update." + config["download_link"]["f_type"])
        params:
            out_dir = join(BASE_FP, "{dataset}", "{dataset}_raw")
        resources:
            time_min=59,
            mem_mb=2048*16
        conda:
            "envs/anndata.yml"
        script:
            "anndata2loom.py" 
        
    rule create_seurat_and_process:
        input:
            counts = join(BASE_FP, "{dataset}", "{dataset}_raw", "X.csv"),
            cell_metadata = join(BASE_FP, "{dataset}", "{dataset}_raw", "obs.csv"),
            gene_metadata = join(BASE_FP, "{dataset}", "{dataset}_raw", "var.csv"),
            h5ad = join(BASE_FP, "{dataset}", "{dataset}_update." + config["download_link"]["f_type"])
        output:
            seurat_temp  =  temp(join(BASE_FP, "{dataset}", "{dataset}_update.h5seurat")), 
            seurat_path = join(BASE_FP, "{dataset}", "seurat_obj.RDS.gz")
        params:
            f_type = config["download_link"]["f_type"],
            use_anndataR = config["use_anndataR"],
            use_seuratdisk = config["use_seuratdisk"],
        conda:  
            "envs/seurat.yml"
        resources:
            time_min=59,
            mem_mb = 96000
        threads: 8
        script:
            "seurat.R"

if config["download_link"]["f_type"] == "GEO":
    rule get_seurat_from_GEO:
        input:
            counts = HTTP.remote(config["download_link"]["GEO"]["counts"]),
            metadata = HTTP.remote(config["download_link"]["GEO"]["metadata"]),
        output:
            seurat_path = join(BASE_FP, "{dataset}", "seurat_obj.RDS.gz")
        params:
            patient_data = config["external_metadata"], 
            geo_data = config["download_link"]["GEO"],
            genes = HTTP.remote(config["download_link"]["GEO"]["genes"])
        conda:
            "envs/seurat.yml"
        resources:
            mem_mb = 140000,
            time_min = "04:59"
        threads: 4
        script:
            "seurat_from_geo.R" 
    
  #  ruleorder: get_seurat_from_GEO > create_seurat_and_process


rule create_cacoa:
    input:
        seurat_path = join(BASE_FP, "{dataset}", "seurat_obj.RDS.gz" )
    output:
        cacoa_obj = join(BASE_FP, "{dataset}", "cao_obj.RDS.gz"),
        umap_coords = join(BASE_FP, "{dataset}", "umap_coords.csv")
    params:
        permute = False
    conda:
        "envs/cacoa.yml"
    threads: 4
    resources:
        mem_mb = 24000,
        time_min = 59
    log:
        "logs/crate_cacoa/{dataset}_create_cacoa.log"
    script:
        "create_cacoa_obj.R"

rule run_cacoa_analysis:
    input:
        cacoa_obj = join(BASE_FP, "{dataset}", "cao_obj.RDS.gz")
    output:
        cacoa_processed =  join(BASE_FP, "{dataset}", "processed_cao.RDS.gz"),
        cao_xlsx = join(BASE_FP, "{dataset}", "report", "{dataset}" + "_de_res.xlsx")
    conda:
        "envs/cacoa.yml"
    threads: 16
    resources:
        mem_mb = 170000,
        time_min = "04:59",
        queue = "long"
    log:
        "logs/run_cacoa/{dataset}_run_cacoa.log"
    script:
        "cacoa_de.R"
    
rule plot_report:
    input:
        cacoa_processed = join(BASE_FP, "{dataset}", "processed_cao.RDS.gz")
    output:
        report_html = join(BASE_FP, "{dataset}", "report", "{dataset}" + "_cacoa.r.ipynb"),
        plot_path = directory(join(BASE_FP, "{dataset}", "report", "{dataset}_plots"))
    params:
        permute = False,
        save_plots = True,
        plot_path = join(BASE_FP, "{dataset}", "report", "{dataset}_plots")
    conda:
        "envs/cacoa.yml"
    log:
        notebook=join(BASE_FP, "{dataset}", "report", "{dataset}" + "_cacoa.r.ipynb")
    resources:
        mem_mb = 192000,
        walltime = "09:59" 
    threads: 4
    notebook:
        "cacoa_report.r.ipynb"

rule export_report:
    input:
        report_notebook = join(BASE_FP, "{dataset}", "report", "{dataset}" + "_cacoa.r.ipynb")
    output:
        report_html =join(BASE_FP, "{dataset}", "report", "{dataset}" + "_cacoa.r.html")
    conda:
        "envs/cacoa.yml"
    log:
        "logs/export_report/{dataset}_export.log"
    shell:
        """
            jupyter nbconvert --to html {input.report_notebook} 
        """
