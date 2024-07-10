# RNA-Seq Analysis walkthrough based on the *npr3/4*, *tga2/5/6*, and *sid2* compared to Col-0 RNA-Seq dataset
## project.sh
[project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh) is a script file built for the GACRC servers and performs the following:
- Checks the quality of raw RNA reads via [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  - [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh) does not perform any sequence cleaning steps (e.g., [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)) as the RNA reads for this instance were cleaned by the sequencing company. We confirm the quality via FastQC.
- Performs a psuedoalignment via [kallisto](https://pachterlab.github.io/kallisto/manual.html) and prepares a directory that can be utilized by [DEG_analysis.R](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/DEG_analysis.R) for differential expression analysis
  - If you want to conduct a traditional alignment analysis, please skip to the [ballgown.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/ballgown.sh) section after reading up to the kallisto section.

The first part of [project.sh](https://github.com/carternewt/RNA_Seq/blob/36b39dfb29bb867d99225d302070ba4fbc6c4406/project.sh) is all of the SLURM headers at the top of the file. 

`#!/bin/bash` This line is necessary everytime and must be present

`#SBATCH --partition=batch` This line is also necessary

`#SBATCH --ntasks=1` This line basically lets you determine if you want to run tasks in your script in parallel. Unless you are doing massive computational analysis this will typically be 1

`#SBATCH --cpus-per-task=12` Depending on how much power you need, you can change the number of CPUs dedicated to your script. Don't be a nuisance to the Sapelo servers and request an unwarranted amount of CPU space. Be realistic, but value your time. 

`#SBATCH --mem=24G` This determines how much memory you need. Typically omic-related analyses are going to need a lot of memory as data has to be stored in the "background" while running programs. 

`#SBATCH --time=3-0:00:00` This determines how long you need the resources you are requesting to run your script. It's always best to overestimate how much time you need then underestimate and have to run your script again. 

`#SBATCH --mail-user=carter.newton@uga.edu` This is to have emails sent to you that update you on how the script is doing

`#SBATCH --mail-type=START,END,FAIL` You can request when you want the servers to update you. In my opinion, it's best to receive notifications when your script finishes and if it fails. 

`#SBATCH --error=/work/lylab/cjn40747/dawei_RNA/project.log.%j.err` This will generate a file in the file path that will contain all information that is outputted by the programs you ran. This is very important to have so you can troubleshoot errors you may encounter in the future. This file will contain error codes that these programs output. Please keep in mind you need to change the file path to a location you are working in. By file path, I am referring to `/work/lylab/cjn40747/dawei_RNA`

`#SBATCH --output=/work/lylab/cjn40747/dawei_RNA/project.log.%j.out` This is similar to the line above. 

You can add other SLURM headers that may be warranted for your analyses and those cna be found [here](https://wiki.gacrc.uga.edu/wiki/Running_Jobs_on_Sapelo2#Header_lines)

---
Next is preparing the data, directories, and programs you'll use. 

`OUT='/work/lylab/cjn40747/dawei_RNA'` This command creates a variable called "OUT" that encodes for the file path found within the quotes. This will make your code more neat and make you less prone to spelling errors as long file paths can be shortened to these variables. You'll want to change the file path to the directory you plan to store most of your files related to the analysis.

`CDNA='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz'` This also creates a variable that now holds a URL to the TAIR10 reference transcriptome I use for this script. Anytime you are able to download your data from a web server, always have your script do it. You'll save space on your own computer and it helps alleviate any errors that may arise with you trying to manually download these large omic files. 

`ml FastQC/0.11.9-Java-11`

`ml kallisto/0.48.0-gompi-2022a`

`ml BWA/0.7.17-GCCcore-11.3.0`

`ml SAMtools/1.16.1-GCC-11.3.0`

`ml BCFtools/1.15.1-GCC-11.3.0`

These are commands that load programs you'll be using within your script. This is also a great way to ensure you have version control/reproducibility. To figure out if GACRC has a program you want, you can use the command `module spider *` and replace "*" with your program name to see if GACRC has uploaded it. 

When it comes to preparing your data, it is likely going to be a fastq file that is gzipped (*.fq.gz). If these can't be downloaded from a web server accessible through your script, upload them to GACRC and move these files into a dedicated directory. 
  - If you have to upload the data to the GACRC servers you can either use the `scp` command through your terminal or use the Globus web server.

In project.sh I uploaded all of the raw fastq RNA reads into the following file path: `/work/lylab/cjn40747/dawei_RNA/all_reads` I also had some reads in another folder, but that's irrelevant to this document.   You should end up with a folder that contains all of your gzipped fastq files with some sort of unique file name. Next, we'll unzip those files 

`gunzip -k $OUT/all_reads/*.fq.gz` 

This is where the variable we created earlier becomes handy. By typing "$OUT" we tell our script to use the value "OUT" encodes and then we can add onto it since I've indicated n teh code above I want it to go to the direcotyr encoded by "OUT" but then into another called "all_reads" and then look for all files that have .fq.gz at the end of their name. the `gunzip` command will unzip all of these .fq.gz files. Prior to running gunzip you should have a folder that contains all your raw RNA reads and will look soemthing like the following: 

![directory of raw reads](https://github.com/carternewt/RNA_Seq/blob/14a05217571fc3d1d3fa1c63cf3e97886123d28d/images/image1.png)

Once the `gunzip` command is run, your folder should be a mix of .fq.gz and .fq files 

![directory of unzipped raw reads](https://github.com/carternewt/RNA_Seq/blob/2c6ed5f94db8b4847cf2cd259a5abcf1ee9eb028/images/image2.png)

Now, we can assess the quality of our reads using FastQC. 

`mkdir -p $OUT/fastqc` First, I created a new directory in which all of my FastQC outputs will be stored.

`fastqc $OUT/all_reads/*.fq -o $OUT/fastqc` Now we run a FastQC analysis on all of our fastq files and use the `-o` argument to tell the FastQC program where we want its generated files stored at. fastQC will generate a lot for us and the most user friendly output is going to be its .html files which can be opened in a web browser. However, viewing all of these html files individually is time-consuming. Instead, we can filter through the summary.txt file that has two main columns. 
1. The first column will either contain the characters "PASS", "WARN", or "FAIL"
2. The second column tells us which test was performed.

We can compile all of these summary.txt files together across all RNA read files and then search for any "FAIL" instances as "WARN" typically is fine to ignore. 

`mkdir -p $OUT/fastqc/all` We're going to create another directory to unzip files into 

`unzip $OUT/fastqc/\*.zip -d $OUT/fastqc/all` FastQC generates zipped folders. So we use the `unzip` command to extract the contents of these folders and then use the `-d` argument to tell the `unzip` command where to extract the files to. 

`find $OUT/fastqc/all -type f -name 'summary.txt' -exec cat {} \; > $OUT/fastqc/all/combined_summary.txt` Now, we'll use the `find` command to search for any file that's called `summary.txt`. Then we use the `-exec` argument to tell the command to "execute" the following command. Thus, we tell our command to `cat` the summary.txt file once it finds it (print the entire content of it) and store it into a new file called "combined_summary.txt". This combined_summary.txt will then be a compiled version of all summary.txt files. 

`grep FAIL $OUT/fastqc/all/combined_summary.txt > $OUT/fastqc/all/fail_summary.txt` Next, we want to find all the "FAIL" instances, so we use `grep.` You can think of `grep` as CTRL + F when you try to search for a phrase in a document or web page. Thus, we use `grep` to find all lines that have "FAIL" in them and then we want to copy those lines into a new file called "fail_summary.txt" so that we can assess what failed the quality assessment more easily. 

If your reads came back with "FAIL" for anything, you may want to check some html files more closely and consider addressing the issues by using Trimmomatic or another cleaning program. 

---

