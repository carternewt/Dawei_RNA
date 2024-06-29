#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=3-0:00:00
#SBATCH --mail-user=carter.newton@uga.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --error=/work/lylab/cjn40747/dawei_RNA/project.log.%j.err
#SBATCH --output=/work/lylab/cjn40747/dawei_RNA/project.log.%j.out

OUT='/work/lylab/cjn40747/dawei_RNA'
CDNA='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz'

ml FastQC/0.11.9-Java-11
ml kallisto/0.48.0-gompi-2022a
ml BWA/0.7.17-GCCcore-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0
ml BCFtools/1.15.1-GCC-11.3.0

gunzip -k $OUT/all_reads/*.fq.gz
gunzip -k $OUT/extra_reads/*.fq.gz
mv $OUT/extra_reads/*.fq $OUT/all_reads
#mkdir -p $OUT/fastqc
#fastqc $OUT/all_reads/*.fq -o $OUT/fastqc
#mkdir -p $OUT/fastqc/all
#unzip $OUT/fastqc/\*.zip -d $OUT/fastqc/all
#find $OUT/fastqc/all -type f -name 'summary.txt' -exec cat {} \; > $OUT/fastqc/all/combined_summary.txt
#grep FAIL $OUT/fastqc/all/combined_summary.txt > $OUT/fastqc/all/fail_summary.txt

curl -s $CDNA | gunzip -c > $OUT/TAIR10.fa
kallisto index -i $OUT/TAIR10.idx $OUT/TAIR10.fa
mkdir -p $OUT/kallisto
for file_1 in $OUT/all_reads/*_1.fq; do
	prefix="${file_1%_1.fq}"
	file_2="${prefix}_2.fq"
	out_dir="$OUT/kallisto/$(basename "$prefix")"
	name=$(basename "$prefix")
	mkdir -p "$out_dir"
	kallisto quant --pseudobam -i $OUT/TAIR10.idx -o $out_dir -b 100 -t 12 "$file_1" "$file_2"
	samtools sort --threads 12 $out_dir/pseudoalignments.bam -o $out_dir/$name.sorted.bam
	samtools index -@ 12 $out_dir/$name.sorted.bam -o $out_dir/$name.sorted.bam.bai
done

mkdir -p $OUT/counts
find $OUT/kallisto -name 'abundance.tsv' -type f | while read -r file; do
	dir=$(dirname "$file")
	dir_name=$(basename "$dir")
	cut -f 1,5 "$file" | awk -v dir="$dir_name" 'BEGIN{OFS="\t"} NR==1 {$2=dir} 1' > "$OUT/counts/$dir_name.tsv"
done
