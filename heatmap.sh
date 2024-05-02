#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --time=8:00:00
#SBATCH --mail-user=carter.newton@uga.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --error=/work/lylab/cjn40747/dawei_RNA/counts/log.%j.err
#SBATCH --output=/work/lylab/cjn40747/dawei_RNA/counts/log.%j.out

source activate heatmap
R --no-save < /home/cjn40747/Dawei_RNA/heatmap_workflow.R
