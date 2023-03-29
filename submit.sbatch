#!/bin/bash -l

#SBATCH -Jtest_pd
#SBATCH --account=hive-jw254
#SBATCH -t30:00:00
#SBATCH -phive
#SBATCH -ooutput.out
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8
##SBATCH -N1 -n8

module load anaconda3/2022.05
conda activate /storage/home/hhive1/jbrandt35/.conda/envs/GalaxiesResearch
conda env config vars set SPS_HOME=/storage/home/hhive1/jbrandt35/data/GalaxiesResearch/Packages/fsps
conda activate /storage/home/hhive1/jbrandt35/.conda/envs/GalaxiesResearch
srun python /storage/home/hhive1/jbrandt35/data/GalaxiesResearch/Packages/powderday/pd_front_end.py /storage/home/hhive1/jbrandt35/data/GalaxiesResearch/Practice parameters_master parameters_model