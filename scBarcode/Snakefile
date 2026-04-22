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
        python3 scripts/handlingTables.py -i "{input.tsv}" -o "{params.outdir}"
        """


def discovered_genotypes(wildcards):
    ckpt = checkpoints.make_barcodes.get(**wildcards)
    barcode_dir = os.path.dirname(ckpt.output.done)
    gts = glob_wildcards(os.path.join(barcode_dir, "{genotype}_barcodes.txt")).genotype
    if not gts:
        raise ValueError(f"No barcode files found in: {barcode_dir}")
    return gts


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
        set -euo pipefail

        mkdir -p "$(dirname "{output.out_bam}")"
        mkdir -p .tmp_subset
        export TMPDIR=.tmp_subset

        subset-bam \
          --bam "{input.bamfile}" \
          --cell-barcodes "{input.barcode_file}" \
          --out-bam "{output.out_bam}" \
          --cores {threads}
        """

rule remove_chrY:
    input:
        bam=f"{config['outdir']}/subset_bams/{{genotype}}.bam"
    output:
        bam=f"{config['outdir']}/subset_bams_noY/{{genotype}}.bam"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120
    singularity:
        "docker://quay.io/biocontainers/samtools:1.22--h96c455f_0"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname "{output.bam}")"

        if [[ "{wildcards.genotype}" == *KOLF* || "{wildcards.genotype}" == *CTL08* ]]; then
            echo "Removing chrY from {wildcards.genotype}"

            samtools index -@ {threads} "{input.bam}"
            refs=$(samtools idxstats "{input.bam}" | cut -f1 | grep -v -E '^chrY$|^Y$' | grep -v '^\*$' | tr '\n' ' ')

            samtools view -@ {threads} -b "{input.bam}" $refs -o "{output.bam}"
            samtools quickcheck "{output.bam}"
        else
            echo "Keeping BAM unchanged for {wildcards.genotype}"
            ln -sf "$(realpath "{input.bam}")" "{output.bam}"
        fi
        """


rule convert_to_fastq:
    input:
        bam=f"{config['outdir']}/subset_bams_noY/{{genotype}}.bam"
    output:
        fastq_dir=directory(f"{config['outdir']}/FASTQs/{{genotype}}"),
        done=touch(f"{config['outdir']}/FASTQs/{{genotype}}/.done")
    params:
        outroot=f"{config['outdir']}/FASTQs"
    singularity:
        "docker://quay.io/biocontainers/10x_bamtofastq:1.4.1"
    threads: 4
    resources:
        runtime=60,
        mem_mb=16000
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outroot}"
        rm -rf "{output.fastq_dir}"

        bamtofastq "{input.bam}" "{output.fastq_dir}"

        test -d "{output.fastq_dir}"
        touch "{output.done}"
        """


rule merge_and_rename_fastqs:
    input:
        fastq_dir=f"{config['outdir']}/FASTQs/{{genotype}}"
    output:
        r1=f"{config['outdir']}/FASTQs/{{genotype}}/{{genotype}}_R1.fastq.gz",
        r2=f"{config['outdir']}/FASTQs/{{genotype}}/{{genotype}}_R2.fastq.gz"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60
    shell:
        r"""
        set -euo pipefail

        r1_files=$(find "{input.fastq_dir}" -type f -name 'bamtofastq_*_R1_*.fastq.gz' | sort -V)
        r2_files=$(find "{input.fastq_dir}" -type f -name 'bamtofastq_*_R2_*.fastq.gz' | sort -V)

        [ -n "$r1_files" ] || { echo "No R1 FASTQs found under {input.fastq_dir}" >&2; exit 1; }
        [ -n "$r2_files" ] || { echo "No R2 FASTQs found under {input.fastq_dir}" >&2; exit 1; }

        cat $r1_files > "{output.r1}"
        cat $r2_files > "{output.r2}"
        """