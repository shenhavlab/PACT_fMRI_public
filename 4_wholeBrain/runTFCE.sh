#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --mem=64GB
#SBATCH --constraint=cascade
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hritz@brown.edu
#SBATCH --account=carney-ashenhav-condo  
#SBATCH -J run-tfce
#SBATCH --output logs/tfce/rdm-tfce-%J.txt

# load module
module load fsl/6.0.3
module load matlab/R2019a


# display info
echo; echo; echo; echo;
echo 'dirPALM: ' 	${dirPALM}
echo 'mask: ' 		${mask}
echo 'dirIN: ' 		${dirIN}
echo 'dirOUT: ' 	${dirOUT}
echo 'nPerm: ' 		${nPerm}
echo; echo; echo; echo;


echo 'started at:'
date
echo; echo; echo; echo;


# matlab-threaded –nodisplay -nodesktop -r "addpath('${dirPALM}'); palm('-i', '${dirIN}', '-m', '${mask}', '-o', '${dirOUT}', '-T', '-accel', 'tail', '-n', ${nPerm}, '-logp')"
matlab-threaded –nodisplay -nodesktop -r "addpath('${dirPALM}'); palm('-i', '${dirIN}', '-m', '${mask}', '-o', '${dirOUT}', '-T', '-n', ${nPerm}, '-logp', '-quiet')"


matlab-threaded –nodisplay -nodesktop -r "flip_nii(${dirOUT}, 1)"


echo; echo; echo; echo;
echo 'finished at:'
date

