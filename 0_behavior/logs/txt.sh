#!/bin/bash
#SBATCH -t 30:00
#SBATCH --mem=128GB
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hritz@brown.edu
#SBATCH --account=carney-ashenhav-condo  
#SBATCH -J senseDyn
#SBATCH --output logs/senseDyn_%J.txt


module load matlab/R2019a



echo 'started at:'
date
echo; echo; echo; echo;


matlab-threaded –nodisplay -nodesktop -r "launch_parpool(12); senseDyn"


echo; echo; echo; echo;
echo 'finished at:'
date
