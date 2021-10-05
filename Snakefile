import pandas as pd
from os.path import join

configfile: "configs/tabularmuris.yaml"

BASE_FP = config["BASE_FP"]
DATASET = config["DATASET"] 

rule all:
    input:
        expand(join(BASE_FP, "{dataset}", "report", "{dataset}" + "_cacoa.r.html"), 
            dataset = DATASET)

rule create_cacoa:
    input:
        seurat_path = join(BASE_FP, "{dataset}", "senis_droplet_lung.rds" )
    output:
        cacoa_obj = join(BASE_FP, "{dataset}", "cao_obj.RDS.gz")
    params:
        permute = False
    conda:
        "envs/cacoa.yml"
    log:
        "logs/crate_cacoa/{dataset}_create_cacoa.log"
    script:
        "create_cacoa.R"

rule run_cacoa_analysis:
    input:
        cacoa_obj = join(BASE_FP, "{dataset}", "cao_obj.RDS.gz")
    output:
        cacoa_processed =  join(BASE_FP, "{dataset}", "processed_cao.RDS.gz")
    conda:
        "envs/cacoa.yml"
    log:
        "logs/run_cacoa/{dataset}_run_cacoa.log"
    script:
        "cacoa_de.R"
    
rule plot_report:
    input:
        cacoa_processed = join(BASE_FP, "{dataset}", "processed_cao.RDS.gz")
    output:
        report_html = join(BASE_FP, "{dataset}", "report", "{dataset}" + "_cacoa.r.ipynb")
 
    params:
        permute = False
    conda:
        "envs/cacoa.yml"
    log:
        "logs/render_report/{dataset}_render_report.log",
        notebook=join(BASE_FP, "{dataset}", "report", "{dataset}" + "_cacoa.r.ipynb")
    resources:
        mem_mb=12000
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
