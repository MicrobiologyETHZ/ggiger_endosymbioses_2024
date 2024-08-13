# To run:
# snakemake -s code/isolate_comparison.py --config norm=False --cluster "DIR=\$(dirname {params.errFile}); mkdir -p \"\${{DIR}}\"; qsub -cwd -S /bin/bash -N snakejob -V -pe smp {threads} -l h_vmem={params.mem} -e {params.errFile} -o {params.outFile}" --jobs 8

import glob,re

samples = [f'GG{i}' for i in range(1,15)] + [f'G4_{str(i).zfill(2)}' for i in range(1,23)]
snps = [f'scratch/{sample}/snippy/{sample}_ref.txt' for sample in samples]
pngs = [f'scratch/{sample}/nucmer/{sample}.png' for sample in samples]
covs = [f'scratch/{sample}/delly/{sample}.cov.gz' for sample in samples]

rule all:
    input:
        snps,
        pngs,
        covs,

rule snippy:
    input:
        ref = "lib/GCF_000198775.1_ASM19877v1_genomic.gbff",
        r1 = "scratch/{sample}/qc/{sample}_R1.fq.gz",
        r2 = "scratch/{sample}/qc/{sample}_R2.fq.gz",
    output:
        summary = "scratch/{sample}/snippy/{sample}_ref.txt",
    params:
        errFile = "scratch/{sample}/snippy/snippy.err.log",
        outFile = "scratch/{sample}/snippy/snippy.out.log",
        outdir = "scratch/{sample}/snippy",
        mem = "4G",
    threads: 32
    shell:
        """
        ml snippy
        snippy --outdir {params.outdir} --ref {input.ref} --R1 {input.r1} --R2 {input.r2} --cpus 32 --prefix {wildcards.sample}_ref --force
        ml purge
        """

rule nucmer:
    input:
        ref = "lib/GCF_000198775.1_ASM19877v1_genomic.fna",
        asm = "scratch/{sample}/assembly/scaffolds.fasta",
    output:
        delta = "scratch/{sample}/nucmer/{sample}.delta",
        png = "scratch/{sample}/nucmer/{sample}.png",
    params:
        errFile = "scratch/{sample}/nucmer/nucmer.err.log",
        outFile = "scratch/{sample}/nucmer/nucmer.out.log",
        mem = "2G",
        prefix = "scratch/{sample}/nucmer/{sample}",
    threads: 32
    shell:
        """
        ml MUMmer
        nucmer {input.ref} {input.asm} -p {params.prefix} --threads 32
        mummerplot --filter --layout --large --png -p {params.prefix} {output.delta}
        ml purge
        """

rule delly:
    input:
        ref = "lib/GCF_000198775.1_ASM19877v1_genomic.fna",
        r1 = "scratch/{sample}/qc/{sample}_R1.fq.gz",
        r2 = "scratch/{sample}/qc/{sample}_R2.fq.gz",
        map = "lib/GCF_000198775.1_ASM19877v1_genomic.map.fa.gz"
    output:
        mark = "scratch/{sample}/delly/{sample}.marked.bam",
        bcf = "scratch/{sample}/delly/{sample}.bcf",
        cov = "scratch/{sample}/delly/{sample}.cov.gz",
    params:
        bam = "scratch/{sample}/delly/{sample}.bam",
        fix = "scratch/{sample}/delly/{sample}.fixed,bam",
        sort = "scratch/{sample}/delly/{sample}.sorted.bam",
        errFile = "scratch/{sample}/delly/delly.err.log",
        outFile = "scratch/{sample}/delly/delly.out.log",
        mem = "4G",
    threads: 32
    shell:
        """
        ml BWA
        ml SAMtools
        bwa mem {input.ref} {input.r1} {input.r2} -t {threads} | samtools view -b > {params.bam}
        samtools fixmate -@ {threads} -m {params.bam} {params.fix}
        samtools sort -@ {threads} {params.fix} > {params.sort}
        samtools markdup -@ {threads} {params.sort} {output.mark}
        samtools index {output.mark}
        rm {params.bam} {params.fix} {params.sort}
        ml delly
        delly call -g {input.ref} {output.mark} > {output.bcf}
        delly cnv -g {input.ref} -m {input.map} {output.mark} -c {output.cov}
        """

