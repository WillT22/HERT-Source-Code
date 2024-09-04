#!/bin/bash

#SBATCH --job-name=HERT_Test
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --partition=general
#SBATCH --time=500:00:00
#SBATCH --mem=20000
#SBATCH --error=error-%j.err
#SBATCH --mail-type=all
#SBATCH --mail-user=shk00082@auburn.edu

module load python/anaconda/3.8.6

HERT.exe Vis.mac
