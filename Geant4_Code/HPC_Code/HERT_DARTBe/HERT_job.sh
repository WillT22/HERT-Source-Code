#!/bin/bash
#SBATCH --job-name=HERT_sim		# job name
#SBATCH --nodes=1			# node(s) required for the job
#SBATCH --ntasks=48			# number of tasks (cores) across all nodes
#SBATCH --partition=mak0037_std		# name of partition to submit job
#SBATCH --output=test-%j.out		# output file where %j is the job ID
#SBATCH --error=test-%j.err		# error file where %j is the job ID
#SBATCH --mail-type=ALL			# will send email for begin, end, and failure
#SBATCH --mail-user=wzt0020@auburn.edu

module load geant4/11.1.2
moudle load cmake

./HERT ../MultipleRun.mac
