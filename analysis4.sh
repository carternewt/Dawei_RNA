#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=1-0:00:00
#SBATCH --mail-user=carter.newton@uga.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --error=/work/lylab/cjn40747/dawei_RNA/analysis4/log.%j.err
#SBATCH --output=/work/lylab/cjn40747/dawei_RNA/analysis4/log.%j.out

source activate analysis4
R --no-save < /home/cjn40747/Dawei_RNA/Analysis4.R
