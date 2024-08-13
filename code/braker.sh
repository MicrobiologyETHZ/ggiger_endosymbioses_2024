#$ -cwd
#$ -S /bin/bash
#$ -N braker
#$ -V
#$ -pe smp 32
#$ -l h_vmem=4G
#$ -e logs/braker.error.log
#$ -o logs/braker.out.log

ml BRAKER

braker.pl --genome=scratch/fungus_2/assembly/assembly.fasta --species=rhizopus_microsporus --bam=data/rnaseq/B1_Aligned.sortedByCoord.out.bam,data/rnaseq/B2_Aligned.sortedByCoord.out.bam --fungus --workingdir=scratch/fungus_2/braker --threads=32 --gff3 --useexisting
