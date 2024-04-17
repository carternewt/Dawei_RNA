#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=3-0:00:00
#SBATCH --mail-user=Dawei.Xu@uga.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --error=/work/lylab/dx99793/Dawei_RNAseq/project.log.%j.err
#SBATCH --output=/work/lylab/dx99793/Dawei_RNAseq/project.log.%j.out

OUT='/work/lylab/dx99793/Dawei_RNAseq'
CDNA='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz'
WGS='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz'

ml FastQC/0.11.9-Java-11
ml SAMtools/1.16.1-GCC-11.3.0
ml Salmon/1.9.0-GCC-11.3.0

gunzip -k $OUT/*.fastq.gz
mkdir -p $OUT/fastqc
fastqc $OUT/*.fastq -o $OUT/fastqc
mkdir -p $OUT/fastqc/all
unzip $OUT/fastqc/\*.zip -d $OUT/fastqc/all
find $OUT/fastqc/all -type f -name 'summary.txt' -exec cat {} \; > $OUT/fastqc/all/combined_summary.txt
grep FAIL $OUT/fastqc/all/combined_summary.txt > $OUT/fastqc/all/fail_summary.txt

curl -s $CDNA > $OUT/TAIR10.fa.gz
gunzip -c $OUT/TAIR10.fa.gz > $OUT/TAIR10.fa
curl -s $WGS > $OUT/TAIR10_WGS.fa.gz
grep '^>' <(gunzip -c $OUT/TAIR10_WGS.fa.gz) | cut -f ' ' -f 1 > $OUT/decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat $OUT/TAIR10.fa.gz $OUT/TAIR10_WGS.fa.gz > $OUT/gentrome.fa.gz
salmon index -t $OUT/gentrome.fa.gz -d $OUT/decoys.txt -p 12 -i $OUT/TAIR10.idx --gencode
mkdir -p $OUT/salmon
for file_1 in $OUT/*_R1_*.fastq; do
	prefix="${file_1%_R1_001.fastq}"
	file_2="${prefix}_R2_001.fastq"
	out_dir="$OUT/salmon/$(basename "$prefix")"
	name=$(basename "$prefix")
	mkdir -p "$out_dir"
	salmon quant -i $OUT/TAIR10.idx -p 12 -l A -1 "$file_1" -2 "$file_2" -o $out_dir/$name -z $out_dir/$name.sam
	samtools view $out_dir/$name.sam -O BAM -o $out_dir/$name.bam
	samtools sort --threads 12 $out_dir/$name.bam -o $out_dir/$name.sorted.bam
	samtools index -@ 12 $out_dir/$name.sorted.bam -o $out_dir/$name.sorted.bam.bai
done
