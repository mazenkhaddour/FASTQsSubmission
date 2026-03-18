import os
from snakemake.io import glob_wildcards
from snakemake.io import directory

configfile: "config.yaml"
def final_outputs(wildcards):
    gts = discovered_genotypes(wildcards)
    return expand(f"{config['outdir']}/FASTQs/{{genotype}}/{{genotype}}_R1.fastq.gz", genotype=gts) + \
           expand(f"{config['outdir']}/FASTQs/{{genotype}}/{{genotype}}_R2.fastq.gz", genotype=gts)


rule all:
    input:
        final_outputs
    default_target: True
wildcard_constraints:
    genotype="[^/]+"

# 1) Checkpoint: create barcode files
checkpoint make_barcodes:
    input:
        tsv=config["cell_ids"]
    output:
        done=touch(f"{config['outdir']}/barcodes/.done")
    params:
        outdir=config["outdir"]
    singularity:
        "docker://mazenkhaddour/demux-script:latest"
    shell:
        r"""
        python3 scripts/handlingTables.py -i {input.tsv} -o {params.outdir}
        """

# 2) Helper: list genotypes after checkpoint runs
def discovered_genotypes(wildcards):
    ckpt = checkpoints.make_barcodes.get(**wildcards)
    barcode_dir = os.path.dirname(ckpt.output.done)
    gts = glob_wildcards(os.path.join(barcode_dir, "{genotype}_barcodes.txt")).genotype
    if not gts:
        raise ValueError(f"No barcode files found in: {barcode_dir}")
    return gts



# 4) Split BAM per genotype
rule split_bam:
    input:
        bamfile=config["bamfile"],
        barcode_file=f"{config['outdir']}/barcodes/{{genotype}}_barcodes.txt",
        done=lambda wc: checkpoints.make_barcodes.get(**wc).output.done
    output:
        out_bam=f"{config['outdir']}/subset_bams/{{genotype}}.bam"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120
    singularity:
        "docker://johnyaku/subset-bam:1.1.0"
    shell:
        r"""
        mkdir -p $(dirname {output.out_bam})
        export TMPDIR=/tmp
        subset-bam \
          --bam {input.bamfile} \
          --cell-barcodes {input.barcode_file} \
          --out-bam {output.out_bam} \
          --cores {threads}
        """
rule convert_to_fastq:
    input:
        bam=f"{config['outdir']}/subset_bams/{{genotype}}.bam"
    output:
        fastq_dir=directory(f"{config['outdir']}/FASTQs/{{genotype}}"),
        done=touch(f"{config['outdir']}/FASTQs/{{genotype}}/.done")
    params:
        outroot=f"{config['outdir']}/FASTQs"
    singularity:
        "docker://quay.io/biocontainers/10x_bamtofastq:1.4.1"
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outroot}

        tmpdir="{params.outroot}/.tmp_bamtofastq_{wildcards.genotype}"
        rm -rf "$tmpdir" "{output.fastq_dir}"

        # IMPORTANT: do NOT mkdir "$tmpdir" — bamtofastq must create it
        bamtofastq "{input.bam}" "$tmpdir"

        # bamtofastq creates the directory we asked for; use it as the produced dir
        if [ ! -d "$tmpdir" ]; then
          echo "ERROR: bamtofastq did not create $tmpdir" >&2
          ls -la {params.outroot} >&2 || true
          exit 1
        fi

        mv "$tmpdir" "{output.fastq_dir}"
        touch "{output.done}"
        """
rule merge_lanes:

    input:
        fastq_dir=f"{config['outdir']}/FASTQs/{{genotype}}"
    output:
        r1=f"{config['outdir']}/FASTQs/{{genotype}}/merged_R1.fastq.gz",
        r2=f"{config['outdir']}/FASTQs/{{genotype}}/merged_R2.fastq.gz"
    shell:
        r"""
        set -euo pipefail

        sample_dir="$(find {input.fastq_dir} -mindepth 1 -maxdepth 3 -type d -name 'Sample_*' | head -n 1)"

        r1s=$(ls "$sample_dir"/bamtofastq_*_R1_001.fastq.gz | sort -V)
        r2s=$(ls "$sample_dir"/bamtofastq_*_R2_001.fastq.gz | sort -V)

        cat $r1s > {output.r1}
        cat $r2s > {output.r2}
        """
rule rename_fastqs:
    input:
        r1=f"{config['outdir']}/FASTQs/{{genotype}}/merged_R1.fastq.gz",
        r2=f"{config['outdir']}/FASTQs/{{genotype}}/merged_R2.fastq.gz",
        fastq_dir=f"{config['outdir']}/FASTQs/{{genotype}}"
    output:
        r1=f"{config['outdir']}/FASTQs/{{genotype}}/{{genotype}}_R1.fastq.gz",
        r2=f"{config['outdir']}/FASTQs/{{genotype}}/{{genotype}}_R2.fastq.gz"
    shell:
        r"""
        set -euo pipefail

        mv {input.r1} {output.r1}
        mv {input.r2} {output.r2}
        """
