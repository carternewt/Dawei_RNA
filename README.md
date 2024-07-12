# RNA-Seq Analysis walkthrough based on the *npr3/4*, *tga2/5/6*, and *sid2* compared to Col-0 RNA-Seq dataset
## In Brief
- [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh)
  - Performs a FastQC quality assessment and pseudoalignment with kallisto
- [ballgown.sh](https://github.com/carternewt/RNA_Seq/blob/da0cf018c4c6b166ec08f87547aa9d53ffe86da6/ballgown.sh)
  - Runs the HISAT2-StringTie-ballgown pipeline, a traditional alignment approach
- [DEG_analysis.R](https://github.com/carternewt/RNA_Seq/blob/c437646f2f9a5ad94dc34b93b89511dea5bd77cc/DEG_analysis.R)
  - An R script that performs a differential expression analysis through edgeR using the output files from [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh)
- [Analysis4.R](https://github.com/carternewt/RNA_Seq/blob/bb6a2e8fda857e53349267003fb3abe3f78c2da8/Analysis4.R) and [analysis4.sh](https://github.com/carternewt/RNA_Seq/blob/bb6a2e8fda857e53349267003fb3abe3f78c2da8/analysis4.sh)
  - Analysis4.R is an R script that is set to be run on the GACRC servers due to it encoding a computationally demanding differential expression analysis. The analysis4.sh script is needed to run Analysis4.R on the GACRC servers.
## project.sh
[project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh) is a script file built for the GACRC servers and performs the following:
- Checks the quality of raw RNA reads via [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  - [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh) does not perform any sequence cleaning steps (e.g., [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)) as the RNA reads for this instance were cleaned by the sequencing company. We confirm the quality via FastQC.
- Performs a psuedoalignment via [kallisto](https://pachterlab.github.io/kallisto/manual.html) and prepares a directory that can be utilized by [DEG_analysis.R](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/DEG_analysis.R) for differential expression analysis
  - If you want to conduct a traditional alignment analysis, please skip to the [ballgown.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/ballgown.sh) section.

The first part of [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh) is all of the SLURM headers at the top of the file. 

`#!/bin/bash` This line is necessary every time and must be present

`#SBATCH --partition=batch` This line is also necessary

`#SBATCH --ntasks=1` This line basically lets you determine if you want to run tasks in your script in parallel. Unless you are doing massive computational analysis, this will typically be 1

`#SBATCH --cpus-per-task=12` Depending on how much power you need, you can change the number of CPUs dedicated to your script. Don't be a nuisance to the Sapelo servers and request an unwarranted amount of CPU space. Be realistic, but value your time. 

`#SBATCH --mem=24G` This determines how much memory you need. Typically omic-related analyses are going to need a lot of memory as data has to be stored in the "background" while running programs. 

`#SBATCH --time=3-0:00:00` This determines how long you need the resources you are requesting to run your script. It's always best to overestimate how much time you need then underestimate and have to run your script again. 

`#SBATCH --mail-user=carter.newton@uga.edu` This is to have emails sent to you that update you on how the script is doing

`#SBATCH --mail-type=START,END,FAIL` You can request when you want the servers to update you. In my opinion, it's best to receive notifications when your script finishes and if it fails. 

`#SBATCH --error=/work/lylab/cjn40747/dawei_RNA/project.log.%j.err` This will generate a file in the file path that will contain all information that is outputted by the programs you ran. This is very important to have so you can troubleshoot errors you may encounter in the future. This file will contain error codes that these programs output. Please keep in mind you need to change the file path to a location you are working in. By file path, I am referring to `/work/lylab/cjn40747/dawei_RNA`

`#SBATCH --output=/work/lylab/cjn40747/dawei_RNA/project.log.%j.out` This is similar to the line above. 

You can add other SLURM headers that may be warranted for your analyses and those can be found [here](https://wiki.gacrc.uga.edu/wiki/Running_Jobs_on_Sapelo2#Header_lines)

---
Next is preparing the data, directories, and programs you'll use. 

`OUT='/work/lylab/cjn40747/dawei_RNA'` This command creates a variable called "OUT" that encodes for the file path found within the quotes. This will make your code more neat and make you less prone to spelling errors as long file paths can be shortened to these variables. You'll want to change the file path to the directory you plan to store most of your files related to the analysis.

`CDNA='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz'` This also creates a variable that now holds a URL to the TAIR10 reference transcriptome I use for this script. Anytime you are able to download your data from a web server, always have your script do it. You'll save space on your own computer and it helps alleviate any errors that may arise with you trying to manually download these large omic files. 

```
ml FastQC/0.11.9-Java-11
ml kallisto/0.48.0-gompi-2022a
ml BWA/0.7.17-GCCcore-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0
ml BCFtools/1.15.1-GCC-11.3.0
```

These are commands that load programs you'll be using within your script. This is also a great way to ensure you have version control/reproducibility. To figure out if GACRC has a program you want, you can use the command `module spider *` and replace "*" with your program name to see if GACRC has uploaded it. 

When it comes to preparing your data, it is likely going to be a FASTQ file that is gzipped (*.fq.gz). If these can't be downloaded from a web server accessible through your script, upload them to GACRC and move these files into a dedicated directory. 
  - If you have to upload the data to the GACRC servers you can either use the `scp` command through your terminal or use the Globus web server.

In project.sh I uploaded all of the raw FASTQ RNA reads into the following file path: `/work/lylab/cjn40747/dawei_RNA/all_reads` I also had some reads in another folder, but that's irrelevant to this document.   You should end up with a folder that contains all of your gzipped FASTQ files with some sort of unique file name. Next, we'll unzip those files 

`gunzip -k $OUT/all_reads/*.fq.gz` 

This is where the variable we created earlier becomes handy. By typing "$OUT," we tell our script to use the value "OUT" encodes, and then we can add onto "OUT." For example, I've indicated in the code above that I want to go to the directory encoded by "OUT" but then into a subdirectory called "all_reads" and then look for all files with .fq.gz at the end of their name. The `gunzip` command will unzip all of these .fq.gz files. Prior to running gunzip, you should have a folder that contains all your raw RNA reads and will look something like the following: 

![directory of raw reads](https://github.com/carternewt/RNA_Seq/blob/14a05217571fc3d1d3fa1c63cf3e97886123d28d/images/image1.png)

Once the `gunzip` command is run, your folder should be a mix of .fq.gz and .fq files 

![directory of unzipped raw reads](https://github.com/carternewt/RNA_Seq/blob/2c6ed5f94db8b4847cf2cd259a5abcf1ee9eb028/images/image2.png)

Now, we can assess the quality of our reads using FastQC. 

`mkdir -p $OUT/fastqc` First, I created a new directory in which all of my FastQC outputs will be stored.

`fastqc $OUT/all_reads/*.fq -o $OUT/fastqc` Now we run a FastQC analysis on all of our FASTQ files and use the `-o` argument to tell the FastQC program where we want its generated files stored at. FastQC will generate a lot for us, and the most user-friendly output is going to be its .html files, which can be opened in a web browser. However, viewing all of these HTML files individually is time-consuming. Instead, we can filter through the summary.txt file that has two main columns. 
1. The first column will either contain the characters "PASS", "WARN", or "FAIL"
2. The second column tells us which test was performed.

We can compile all of these summary.txt files together across all RNA read files and then search for any "FAIL" instances. 

`mkdir -p $OUT/fastqc/all` We're going to create another directory to unzip files into 

`unzip $OUT/fastqc/\*.zip -d $OUT/fastqc/all` FastQC generates zipped folders. So we use the `unzip` command to extract the contents of these folders and then use the `-d` argument to tell the `unzip` command where to extract the files to. 

`find $OUT/fastqc/all -type f -name 'summary.txt' -exec cat {} \; > $OUT/fastqc/all/combined_summary.txt` Now, we'll use the `find` command to search for any file that's called `summary.txt.` Then, we use the `-exec` argument to tell the command to "execute" the following command. Thus, we tell our command to `cat` the summary.txt file once it finds it (print the entire content of it) and store it in a new file called "combined_summary.txt". This combined_summary.txt will then be a compiled version of all summary.txt files. 

`grep FAIL $OUT/fastqc/all/combined_summary.txt > $OUT/fastqc/all/fail_summary.txt` Next, we want to find all the "FAIL" instances, so we use `grep.` You can think of `grep` as CTRL + F when you try to search for a phrase in a document or web page. Thus, we use `grep` to find all lines that have "FAIL" in them, and then we want to copy those lines into a new file called "fail_summary.txt" so that we can assess what failed the quality assessment more easily. You can then go to the directory where you saved your fail_summary.txt file to and use the `cat` command to view all "FAIL" instances. 

![fail_summary.txt content](https://github.com/carternewt/RNA_Seq/blob/da0cf018c4c6b166ec08f87547aa9d53ffe86da6/images/image3.png)

If your reads came back with "FAIL" for anything, you might want to check some HTML files more closely and consider addressing the issues by using Trimmomatic or another cleaning program. 

---
Next is getting files ready for kallisto. 

`curl -s $CDNA | gunzip -c > $OUT/TAIR10.fa` For kallisto we need a reference that our reads can be mapped to. kallisto utilizes a transcriptome (cDNA) library for its reference. The variable created earlier in this script, "CDNA," contains a link to the Ensembl database where the FASTA file for the TAIR10 cDNA library is. By using `curl`, we can obtain data from URLs such as the one encoded by "CDNA." Once we have downloaded the TAIR10 FASTA file, we need to unzip it as we did for the FASTQ files earlier, as this file is in .fa.gz format.

`kallisto index -i $OUT/TAIR10.idx $OUT/TAIR10.fa` Now, we need to let Kkllisto read through the FASTA file and index it so that the program can utilize the FASTA file for referencing. We use the `-i` argument here to tell `kallisto index` what we want to name the index file and where to store it. 

`mkdir -p $OUT/kallisto` We're also going to create a directory where we want our files to be stored during the kallisto analysis. 

```
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
```
This for loop is the most important part of this script. Since we have multiple RNA reads to run through, we want to make sure the script can run through all of them instead of having to run each set of RNA reads separately. Thus, making sure your files are organized is **VERY** important. We're going to break down lines of code into segments as needed. 

**NOTE:** I'll use Col-0-1 as an example to help conceptualize codes for the following lines. 

For this for loop to function, we want it to look in a folder where we've stored all of our RNA reads at that we want to run through kallisto. `for file_1 in $OUT/all_reads/*_1.fq` here we are creating a new variable (similar to "OUT" and "CDNA" called "file_1". "file_1" is going to encode for a file that ends in _1.fq that is found in the all_reads directory. Thus, we are having our for loop find one of the paired reads. Then, we end the for loop with `; do`, which tells it to now carry out all the following lines of code once it finds a file that fits the criteria. Make sure you indent all following lines of code. 

The next couple of lines consist of creating variables and directories so that we can keep our kallisto outputs organized. `prefix="${file_1%_1.fq}"` is creating a variable called "prefix" that is going to be equal to the name of a file up until _1.fq. For example, if the file the for loop is working with is "Col-0-1_1.fq", then prefix will be equal to "Col-0-1". 

Next, we need to find the other paired read for the sample we are working with as we only know where one is. `file_2="${prefix}_2.fq"` creates a new variable called "file_2" and finds the other paired read file by utilizing the "prefix" variable we created earlier. For example, assuming "prefix" is still equal to "Col-0-1" that means "file_2 will be equal to "Col-0-1_2.fq" 

Now, we want to create a directory specifically for "Col-0-1" that we can output or kallisto files into. `out_dir="$OUT/kallisto/$(basename "$prefix")"` To do this, we'll create another variable that's going to be equal to the file path for this specific Col-0-1 directory. The `$(basename "$prefix")` code is a bit to understand. We use parentheses here because `basename` is a command that prints out variables. We want to run the `basename` command on our variable "prefix" before the other lines in our code, so we must put it in parentheses. We also always include "$" before variables so that the code actually uses the information it encodes for. Thus, the `$(basename "$prefix")` command will output "Col-0-1". Put that all together, and this code will create a new variable that encodes a file path to "$OUT/kallisto/Col-0-1".

Now we need to create this new file path we just generated. `mkdir -p "$out_dir" ` will then create a new directory using the value encoded by "out_dir" which we just created previously. Thus, we will create a new directory at "$OUT/kallisto/Col-0-1". 

Everything is now prepped to run kallisto since we have a unique directory to store our files in. The kallisto command is fairly straightforward. We're going to call on `kallisto quant` to run the quantification algorithm. The following are additional arguments we need to supply to `kallisto quant` for it to function. 
- `--pseudobam` is an optional argument that only exists in kallisto versions older than version 0.48.0. This argument will generate a BAM file of the pseudoalignments, which can be used to visualize alignments through an optional downstream analysis if wanted.
- `-i $OUT/TAIR10.idx` is an argument that tells kallisto where it can find the index file we generated with `kallisto index` previously.
- `-o $out_dir` is an argument to inform kallisto where to output the files it generates.
- `-b 100` is an argument that tells kallisto the number of bootstrap samples. Bootstrapping is important to increase confidence and reliability and identify biases in the quantification step.
- `-t 12` is an argument that indicates the number of threads Kallisto can utilize. Make this value equal to the number of CPUs you requested for your script, which can be found in the SLURM headers.
- `"$file_1" "$file_2"` once we've added all of our additional arguments, we need to supply the paired reads we want kallisto to utilize. Remember, our for loop generated the "file_1" and "file_2" variables that encode the names of our FASTQ files.

```
samtools sort --threads 12 $out_dir/pseudoalignments.bam -o $out_dir/$name.sorted.bam
samtools index -@ 12 $out_dir/$name.sorted.bam -o $out_dir/$name.sorted.bam.bai
```
These next two lines of code aren't necessary unless you want to visualize how the reads are aligning to the reference. In brief, these lines will take the BAM file generated by kallisto (if you use the --pseudobam argument) and generate an indexed BAM file that can be used for visualisation.

**NOTE:** Make sure you end your for loop with `done` and that it isn't indented. 

---
With the kallisto quantification done, we simply need to organize our output so that it's in a format that can be read for subsequent analyses. In this case, the [DEG_analysis.R](https://github.com/carternewt/RNA_Seq/blob/c437646f2f9a5ad94dc34b93b89511dea5bd77cc/DEG_analysis.R) script was used to analyze the kallisto outputs. To read in data for this R script, we need the "abundance.h5" files generated by kallisto. We also can't rename these "abundance.h5" files, as this will give us issues when importing our data into the R script. Instead, we just need to make sure the abundance.h5 files are in a directory that's named in accordance with the sample it represents, and all of these directories are housed in another directory. The final directory we are aiming to make will be structured like this:

```
$OUT/h5_files/
  Col-0-1/
    abundance.h5
  Col-0-2/
    abundance.h5
  Col-0-3/
    abundance.h5
...
```

To achieve this, the following code will be employed
```
mkdir -p $OUT/h5_files
find $OUT/kallisto -name 'abundance.h5' -type f | while read -r file; do
        dir=$(dirname "$file")
        out_dir="$OUT/h5_files/$(basename "$dir")"
        mkdir -p $out_dir
        cp "$file" $out_dir
done
```

First, we need to create the directory that's going to house all of our subdirectories. `mkdir -p $OUT/h5_files` is where we will create all of our subdirectories that will hold respective abundance.h5 files. 

Next, we are going to create another looping code. Instead of using `for` to generate a looping command, we are going to use `while`. It's essentially the same idea as both of these commands will continue to run until the conditions we give them have run out. 

To start our `while` loop we are going to provide the condition for it. Our condition is going to be finding all of the abundance.h5 files. We'll be using the `find` command we utilized earlier in the script to do so. `find $OUT/kallisto -name 'abundance.h5' -type f` here, we are telling the find command to search for files within the kallisto directory that have file name "abundance.h5". Once it has found all those files we want to input all the found files into the `while` command which is why we have a `|` between the commands. The `|` acts as a pipe and tells your script to take the output from the code on the left side of the pipe and use it as input for the code on the right side of the pipe. Then our "looping" command `while read -r file; do` is stating to read the "file" (these are the abundance.h5 files) one at a time. With these files we essentially want to move them into a directory to align with our formatting. 

```
dir=$(dirname "$file")
out_dir="$OUT/h5_files/$(basename "$dir")"
mkdir -p $out_dir
```

These three lines should be recognizable from the `for` loop earlier, so I will skip over this. The only new command is the following one, which is `cp "$file" $out_dir` where we are using `cp` to copy the abundance.h5 file into its new directory.

That is it for project.sh! You should be able to navigate in Sapelo to the directories where you've saved all of these files. You'll need to end up copying certain files or folders (such as the h5_files directory we just created) onto your local computer for downstream analyses. As a reminder, this can be achieved with `scp` or Globus. 

## ballgown.sh
[ballgown.sh](https://github.com/carternewt/RNA_Seq/blob/da0cf018c4c6b166ec08f87547aa9d53ffe86da6/ballgown.sh) is a script file built for the GACRC servers and performs the following: 
- Employs [HISAT2](https://daehwankimlab.github.io/hisat2/) to align reads to a reference genome
- Pipes HISAT2 alignments into [StringTie](https://ccb.jhu.edu/software/stringtie/) to assemble potential transcripts
  - Generates estimated transcript abundance files (.ctab) that can be used by [ballgown](https://github.com/alyssafrazee/ballgown) to read into R
 
The first part of [ballgown.sh](https://github.com/carternewt/RNA_Seq/blob/da0cf018c4c6b166ec08f87547aa9d53ffe86da6/ballgown.sh) consists of SLURM headers similar to [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh). For an explanation on these lines, please refer to the start of the project.sh section. 

The first 11 lines of [ballgown.sh](https://github.com/carternewt/RNA_Seq/blob/da0cf018c4c6b166ec08f87547aa9d53ffe86da6/ballgown.sh) are similar to [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh) in regards to the commands used. 

```
OUT='/work/lylab/cjn40747/dawei_RNA'
REF='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_rm.toplevel.fa.gz'
GFF3='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gff3.gz'
```
HISAT2 requires a genome FASTA file. Thus, the variable "REF" encodes for an *Arabidopsis thaliana* TAIR10 reference genome. StringTie requires a general feature format (GFF or GFF3) file that contains annotations for genomic features. Thus, the variable "GFF#" encodes for the *Arabidopsis thaliana* TAIR10 GFF3 file. 

```
ml HISAT2/3n-20201216-gompi-2022a
ml SAMtools/0.1.20-GCC-11.2.0
ml StringTie/2.2.1-GCC-11.3.0 
ml gffread/0.12.7-GCCcore-11.3.0
```

Here, we load in the programs that will be utilized in this script. 

---
Next, we need to prep our data. 
```
curl -s $REF | gunzip -c > $OUT/TAIR10_DNA.fa
hisat2-build -f -p 12 $OUT/TAIR10_DNA.fa $OUT/TAIR10_DNA_idx
curl -s $GFF3 | gunzip -c > $OUT/TAIR10_DNA.gff
gffread $OUT/TAIR10_DNA.gff -T -o $OUT/TAIR10_DNA.gtf
```
Similar to [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh), to obtain the files from the Ensembl database, we'll utilize the `curl` command piped into `gunzip` to unzip the files. 

To prepare the genomic sequence FASTA file for HISAT2, the program needs an indexed version it can read, similar to what was performed for kallisto. `hisat2-build -f -p 12 $OUT/TAIR10_DNA.fa $OUT/TAIR10_DNA_idx` this command will generate the indexed reference genome for us to be called upon later in the script. The `-f` argument tells `hisat2-build` that the input file we are giving it is in FASTA format. `-p 12` tells `hisat2-build` to utilize 12 CPUs. Change "12" to however many CPUs you requested in your script which is indicated in the SLURM headers. 

For StringTie to assemble reads into transcripts, it requires a GTF file. Since Ensembl provides a GFF file, we need to convert the GFF file to a GTF. `gffread $OUT/TAIR10_DNA.gff -T -o $OUT/TAIR10_DNA.gtf` performs this conversion for us and doesn't require many arguments. `-o` is the only main argument needed, which is to tell `gffread` what we want to name out GTF file and where to store it. 

---
With our data prepared, it's time to run our reads through the HISAT2-StringTie pipeline. It is also assumed that your RNA reads are all stored in a folder and unzipped, as discussed previously in the [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh) section. 

```
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
```

The first 7 lines are identical to the `for` loop in [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh) that was used to run the kallisto quantifcation. The only difference is that we have created a new directory for our files to be stored in. That said, we'll skip ahead to the new lines of code where we implement commands from HISAT2, SAMtools, and StringTie.

`hisat2 -p 12 --dta -x $OUT/TAIR10_DNA_idx -1 "$file_1" -2 "$file_2" -S $out_dir/$name.sam` the first command we'll use is `hisat2` which will align our reads to the indexed reference genome we generated and output the alignments into a SAM file. `-p 12` tells the `hisat2` command how many CPUs to use. `--dta` is an argument that ensures the reported alignments are tailored for StringTie to utilize. `-x` argument is used to indicate the **basename** of the indexed reference genome. By basename, this means we do not indicate the file type. This is because when we generated the indexed reference genome, it generated multiple .ht2 files, as seen below.

![Indexed reference genome](https://github.com/carternewt/RNA_Seq/blob/d14ceb35302f9c8749c85baf3e69a84ab9e5391a/images/image4.png)

Thus, we only call on the basename to ensure that `hisat2` can find all the indexed files for utilization. `-1 "$file_1" -2 "$file_2"` the `-1` and `-2` arguments have to be followed by the paired reads we are providing `hisat2`. Similar to [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh), we set the variables "file_1" and "file_2" equal to those paired reads. `-S` argument is used to tell `hisat2` what to name the SAM file it will generate and where it should be stored. 

Next, we need to convert the SAM file generated by HISAT2 into a BAM file and sort it. We'll use two different SAMtools commands to achieve this. First, `samtools view -bS $out_dir/$name.sam -o $out_dir/$name.bam` takes our SAM file from HISAT2 and converts it into a BAM file. `samtools view` is the command we'll utilize and requires only a couple of arguments. `-bS` is two arguments combined into one. The "b" component tells `samtools view` to output a BAM file. The "S" component tells `samtools view` that we are providing it a SAM file for input. If you use a more recent version of SAMtools, consult the manual to see if the `-S` argument is still necessary. `-o` argument tells `samtools view` what to name the generated BAM file and where to store it. 

With our BAM file generated, we'll need to sort it so that StringTie can properly utilize it. `samtools sort -@ 12 $out_dir/$name.bam $out_dir/$name.sorted` the code for sorting a BAM file is fairly straightforward. The only argument needed is `-@`, where we tell `samtools sort` how many CPUs it can use. After this, we supply the location of the BAM file, followed by where we want the sorted BAM file to be stored. It's important to follow this ordering of input first then output. 

Now, we can utilize StringTie to assemble potential transcripts. `stringtie -p 12 -G $OUT/TAIR10_DNA.gtf -o $out_dir/$name.gtf -l $name $out_dir/$name.sorted.bam` is also fairly straightforward. `-p` tells `stringtie` how many CPUs it can utilize. `-G` tells `stringtie` where it can find the GTF file to use as an annotation reference. `-o` argument tells `stringtie` what to name the generated GTF file and where to store them. `-l` argument provides `stringtie` the basename to use as a prefix for output transcripts. After we supply all of these arguments, we provide `stringtie` with the location of the sorted BAM file we want it to use for assembling transcripts. 

---
The `for` loop above will have generated all of the files that are needed for downstream analyses. Thus, we need to organize relevant data into a format that downstream programs can readily interpret (e.g., ballgown and edgeR). 




