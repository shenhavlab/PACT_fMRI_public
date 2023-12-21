#!/bin/bash
#SBATCH -t 30:00
#SBATCH --mem=128GB
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hritz@brown.edu
#SBATCH --account=carney-ashenhav-condo  
#SBATCH -J rdm-collect
#SBATCH --output logs/rdm-collect_level1_%J.txt


module load matlab/R2019a



echo 'started at:'
date
echo; echo; echo; echo;


matlab-threaded –nodisplay -nodesktop -r "RDM_collect_level1;"


echo; echo; echo; echo;
echo 'finished at:'
date
