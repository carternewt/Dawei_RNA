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
REF='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_rm.toplevel.fa.gz'
GFF3='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gff3.gz'

ml HISAT2/3n-20201216-gompi-2022a
ml SAMtools/0.1.20-GCC-11.2.0
ml StringTie/2.2.1-GCC-11.3.0 #Feng used StringTie/2.1.1-GCC-8.30 but this version is not on Sapelo anymore
ml gffread/0.12.7-GCCcore-11.3.0

curl -s $REF | gunzip -c > $OUT/TAIR10_DNA.fa
hisat2-build -f -p 12 $OUT/TAIR10_DNA.fa $OUT/TAIR10_DNA_idx
curl -s $GFF3 | gunzip -c > $OUT/TAIR10_DNA.gff
gffread $OUT/TAIR10_DNA.gff -T -o $OUT/TAIR10_DNA.gtf

mkdir -p $OUT/ballgown_version
for file_1 in $OUT/all_reads/*_1.fq; do
	prefix="${file_1%_1.fq}"
	file_2="${prefix}_2.fq"
	out_dir="$OUT/ballgown_version/$(basename "$prefix")"
	name=$(basename "$prefix")
	mkdir -p "$out_dir"
	hisat2 -p 12 --dta -x $OUT/TAIR10_DNA_idx -1 "$file_1" -2 "$file_2" -S $out_dir/$name.sam
	samtools view -bS $out_dir/$name.sam -o $out_dir/$name.bam
	samtools sort -@ 12 $out_dir/$name.bam $out_dir/$name.sorted
	stringtie -p 12 -G $OUT/TAIR10_DNA.gtf -o $out_dir/$name.gtf -l $name $out_dir/$name.sorted.bam
done

find $OUT/ballgown_version -type f -name "*.gtf" > $OUT/ballgown_version/stringtie_merge.txt
stringtie --merge -p 12 -G $OUT/TAIR10_DNA.gtf -o $OUT/ballgown_version/stringtie_merge.gtf $OUT/ballgown_version/stringtie_merge.txt

mkdir -p $OUT/ballgown_version/sortedBAM
find $OUT/ballgown_version -name "*.sorted.bam" -type f | while read -r file; do
	mv "$file" $OUT/ballgown_version/sortedBAM
done

mkdir -p $OUT/ballgown_version/ballgown_input_files
for file in $OUT/ballgown_version/sortedBAM/*.sorted.bam; do
	prefix="${file%.sorted.bam}"
	name=$(basename "$prefix")
	out_dir="$OUT/ballgown_version/ballgown_input_files/$(basename "$prefix")"
	mkdir -p $out_dir
	stringtie -e -B -p 12 -G $OUT/ballgown_version/stringtie_merge.gtf -o $out_dir/$name.gtf "$file"
done
