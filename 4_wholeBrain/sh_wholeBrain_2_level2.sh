#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH --mem=128GB
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hritz@brown.edu
#SBATCH --account=carney-ashenhav-condo  
#SBATCH -J rdm-lvl2
#SBATCH --output logs/rdm-lvl2-%J.txt




#--------- Variables ---------

# root directory
root_dir="/users/hritz/data/mri-data/RDM2"
spm_dir="/users/hritz/data/mri-data/analysistools/spm12"


# name
name="feature" 		# 

# FWHM
fwhm=8				#

# confound
confound="movt" 	#

# AR
cvi="wls" 			# 




module load matlab/R2019a
module load spm/spm12


echo; echo; echo; echo;
echo "========== rdm lvl 2 =========="
echo "dir: '${root_dir}' | name: '${name}' | fwhm: ${fwhm} | confound: '${confound}' | cvi: '${cvi}'"
echo; echo; echo; echo;

echo 'started at:'
date
echo; echo; echo; echo;

# do level 2 estimate
matlab-threaded –nodisplay -nodesktop -r "RDM_wholeBrain_2_1_estimate('${root_dir}', '${spm_dir}', '${name}', ${fwhm}, '${confound}', '${cvi}')"

echo; echo; echo; echo;
echo 'finished at:'
date





echo; echo; echo; echo;
echo "========== rdm lvl 2 contrast =========="
echo "dir: '${root_dir}' | name: '${name}' | fwhm: ${fwhm} | confound: '${confound}' | cvi: '${cvi}'"
echo; echo; echo; echo;

echo 'started at:'
date
echo; echo; echo; echo;

# do level 2 contrast
matlab-threaded –nodisplay -nodesktop -r "RDM_wholeBrain_2_2_contrast('${root_dir}', '${spm_dir}', '${name}',  ${fwhm}, '${confound}', '${cvi}')"

echo; echo; echo; echo;
echo 'finished at:'
date
