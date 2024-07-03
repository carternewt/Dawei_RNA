#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=3-0:00:00
#SBATCH --mail-user=carter.newton@uga.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --error=/work/lylab/cjn40747/dawei_RNA/ballgown_version/project.log.%j.err
#SBATCH --output=/work/lylab/cjn40747/dawei_RNA/ballgown_version/project.log.%j.out

OUT='/work/lylab/cjn40747/dawei_RNA'
REF='

ml HISAT2/3n-20201216-gompi-2022a
ml SAMtools/0.1.20-GCC-11.2.0
ml StringTie/2.2.1-GCC-11.3.0 #Feng used StringTie/2.1.1-GCC-8.30 but this version is not on Sapelo anymore

mkdir -p $OUT/ballgown_version/hisat
for file_1 in $OUT/all_reads/*_1.fq; do
	prefix="${file_1%_1.fq}"
	file_2="${prefix}_2.fq"
	out_dir="$OUT/ballgown_version/hisat/$(basename "$prefix")"
	name=$(basename "$prefix")
	mkdir -p "$out_dir"
	
