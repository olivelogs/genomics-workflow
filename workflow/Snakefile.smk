rule fastqc:
    input:
        directory("resources/data/)
    output:
        "results/fastqc/{sample}_fastqc.zip"
        "results/fastqc/{sample}_fastqc.html"
    params:
        OUTPUT_DIR="results/fastqc",
        THREADS=4
    shell:
        """
        bash scripts/01_fastqc.sh {INPUT_DIR} {params.OUTPUT_DIR} {params.THREADS}
        """

rule process_radtags:
    input:
        "resources/data/{sample}.fq"
    output:
        "results/clean/{sample}_clean.fq"
    params:
        OUTPUT_DIR="results/radtags",
        BARCODE=""
    shell:
        """
        bash scripts/02_process_radtags.sh {input.fq} {params.OUTPUT_DIR} {params.BARCODE}
        """

# DAG #
SAMPLES = ["sample1","sample2","sample3"]  # or load from config/samples.txt

rule all:
    input:
        expand("results/radtags/{sample}.clean.fq", sample=SAMPLES)