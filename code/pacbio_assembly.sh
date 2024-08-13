#$ -cwd
#$ -S /bin/bash
#$ -N flye
#$ -V
#$ -pe smp 32
#$ -l h_vmem=8G
#$ -e logs/flye.error.log
#$ -o logs/flye.out.log

for sample in 1 2
do
mkdir -p scratch/fungus_$sample/
mkdir -p scratch/fungus_$sample/assembly

ml Flye/2.9.2-GCC-10.2.0

flye --pacbio-hifi data/pacbio/BMK230222-BH052-02P000$sample\.ccs.fastq.gz -o scratch/fungus_$sample/assembly/ -t 32
done
