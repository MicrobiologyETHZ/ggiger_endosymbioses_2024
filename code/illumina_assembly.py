# To run:
# snakemake -s code/isolate_snake.py --config norm=False --cluster "DIR=\$(dirname {params.errFile}); mkdir -p \"\${{DIR}}\"; qsub -cwd -S /bin/bash -N snakejob -V -pe smp {threads} -l h_vmem={params.mem} -e {params.errFile} -o {params.outFile}" --jobs 8

import glob,re

samples = [f'GG{i}' for i in range(1,15)] + [f'R00{str(i).zfill(2)}' for i in range(1,23)] + [f'G4_{str(i).zfill(2)}' for i in range(1,23)]
asms = [f'scratch/{sample}/assembly/{sample}.ass.done' for sample in samples]
motus = [f'scratch/{sample}/motus/{sample}.motus' for sample in samples]

rule all:
    input:
        asms,
        motus,

rule qc:
    input:
        r1 = "data/raw/{sample}_R1.fq.gz",
        r2 = "data/raw/{sample}_R2.fq.gz"
    output:
        qc1 = "scratch/{sample}/qc/{sample}_R1.fq.gz",
        qc2 = "scratch/{sample}/qc/{sample}_R2.fq.gz",
        fail = "scratch/{sample}/qc/{sample}_fail.fq.gz",
        single = "scratch/{sample}/qc/{sample}_single.fq.gz",
        stats = "scratch/{sample}/qc/{sample}_qc.stats",
        adaptMatch = "scratch/{sample}/qc/{sample}_adaptMatch.fq.gz",
        adaptSingle = "scratch/{sample}/qc/{sample}_adaptSingle.fq.gz",
        adaptStats = "scratch/{sample}/qc/{sample}_adapt.stats",
        phixMatch = "scratch/{sample}/qc/{sample}_phixMatch.fq.gz",
        phixSingle = "scratch/{sample}/qc/{sample}_phixSingle.fq.gz",
        phixStats = "scratch/{sample}/qc/{sample}_phix.stats"
    params:
        errFile = "scratch/{sample}/qc/qc.err.log",
        outFile = "scratch/{sample}/qc/qc.out.log",
        mem = "2G"
    threads: 32
    log:
        log = "scratch/{sample}/qc/{sample}.qc.log",
        adaptLog = "scratch/{sample}/qc/{sample}.adapt.log",
        phixLog = "scratch/{sample}/qc/{sample}.phix.log"
    shell:
        """
        ml BBMap
        bbduk.sh -Xmx2G usejni=t threads=14 overwrite=t qin=33 in1={input.r1} in2={input.r2} ref=/nfs/modules/modules/software/BBMap/38.26-foss-2018b/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 out=stdout.fq outm={output.adaptMatch} outs={output.adaptSingle} refstats={output.adaptStats} statscolumns=5 2> {log.adaptLog} | bbduk.sh -Xmx2G usejni=t threads=4 interleaved=true overwrite=t qin=33 in=stdin.fq out=stdout.fq outm={output.phixMatch} outs={output.phixSingle} ref=/nfs/modules/modules/software/BBMap/38.26-foss-2018b/resources/phix174_ill.ref.fa.gz k=31 hdist=1 refstats={output.phixStats} statscolumns=5 2> {log.phixLog} | bbduk.sh -Xmx2G usejni=t threads=14 overwrite=t interleaved=true qin=33 in=stdin.fq fastawrap=10000 out1={output.qc1} out2={output.qc2} interleaved=true outm={output.fail} outs={output.single} minlength=45 maq=20 maxns=0 overwrite=t stats={output.stats} statscolumns=5 trimq=14 qtrim=r 2>{log.log}
        ml purge
        """

rule assembly:
    input:
        qc1 = "scratch/{sample}/qc/{sample}_R1.fq.gz",
        qc2 = "scratch/{sample}/qc/{sample}_R2.fq.gz",
        single = "scratch/{sample}/qc/{sample}_single.fq.gz"
    output:
        done = touch("scratch/{sample}/assembly/{sample}.ass.done"),
    params:
        outdir = "{sample}",
        errFile = "scratch/{sample}/assembly/spades.err.log",
        outFile = "scratch/{sample}/assembly/spades.out.log",
        mem = "4G"
    threads: 32
    log:
        log = "scratch/{sample}/assembly/{sample}.ass.log"
    shell:
        """
        ml SPAdes
        spades.py --isolate --pe1-1 {input.qc1} --pe1-2 {input.qc2} --pe1-s {input.single} -o scratch/{params.outdir}/assembly -t {threads} -m 128
        ml purge
        """

rule motus:
    input:
        qc1 = "scratch/{sample}/qc/{sample}_R1.fq.gz",
        qc2 = "scratch/{sample}/qc/{sample}_R2.fq.gz",
        single = "scratch/{sample}/qc/{sample}_single.fq.gz"
    output:
        results = "scratch/{sample}/motus/{sample}.motus",
    params:
        errFile = "scratch/{sample}/motus/motus.err.log",
        outFile = "scratch/{sample}/motus/motus.out.log",
        mem = "4G",
        name = "{sample}"
    threads: 32
    log:
        log = "scratch/{sample}/motus/{sample}.motus.log"
    shell:
        """
        ml mOTUs
        motus profile -f {input.qc1} -r {input.qc2} -s {input.single} -n {params.name} -o {output.results} -t {threads}
        ml purge
        """ 

#rule filter:
#    input:
#        done = "assembly/{sample}.ass.done",
#        scaffolds = "assembly/scaffolds.fasta"
#    output:
#        scaffolds = "assembly/{sample}_scaffolds_1k.fasta"
#    params:
#        errFile = "filter.err.log",
#        outFile = "filter.out.log",
#        pe = "smp 1",
#        mem = "8G"
#    shell:
#        """
#        ml Perl
#        perl /nfs/home/fieldc/scripts/contig-filter.pl -f {input.scaffolds} -c 1000 > {output.scaffolds}
#        ml purge
#        """

#rule prodigal:
#    input:
#        scaffolds = "assembly/scaffolds.fasta"
#    output:
#        fna = "annotation/{sample}.fna",
#        faa = "annotation/{sample}.faa",
#        gff = "annotation/{sample}.gff",
#        gbk = "annotation/{sample}.gbk"
#    params:
#        errFile = "",
#        outFile = "",
#        pe = "smp 16",
#        mem = "4G"
#    shell:
#        """
#        ml prodigal
#        prodigal -a {output.faa} -d {output.fna} -f gff -o {output.gff} -c -q -m -i {input.scaffolds}
#        prodigal -o {output.gbk} -c -q -m -i {input.scaffolds}
#        ml purge
#        """

#rule emapper:
#    input:
#        faa = "annotation/{sample}.faa"
#    output:
#        hmm = "annotation/{sample}.emapper.hmm_hits",
#        seed = "annotation/{sample}.emapper.seed_orthologs",
#        anno = "annotation/{sample}.emapper.annotations"
#    params:
#        errFile = "annotation/emapper.err.log",
#        outFile = "annotation/emapper.out.log",
#        pe = "smp 16",
#        mem = "4G",
#        prefix = "{sample}"
#    shell:
#        """
#        ml eggnog-mapper/1.0.3-foss-2018b-Python-2.7.15
#        emapper.py -d bact -o annotation/{params.prefix} -i {input.faa} --cpu 16
#        ml purgeq
#        """
