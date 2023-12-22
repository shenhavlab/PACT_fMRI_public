#!/bin/bash
#SBATCH -t 0:30:00
#SBATCH --mem=32GB
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hritz@brown.edu
#SBATCH --account=carney-ashenhav-condo  
#SBATCH -J rdm-eigenvariate_pt
#SBATCH --output logs/rdm-eigenvariate_%A-%a.txt
#SBATCH --array=4-32



#--------- Variables ---------

# root directory
root_dir="/users/hritz/data/mri-data/RDM2"
spm_dir="/users/hritz/data/mri-data/analysistools/spm12"


# name
name="perfCongBlk"      # 

# FWHM
fwhm=0				# 

# confound
confound="movt" 	# 

# AR
cvi="wls" 			#



#----------- Dictionaries for subject specific variables -----
# Dictionary of labels per subject             
declare -A labels=(\
                    [4]="8004" \
                    [5]="8005" \
                    [6]="8006" \
                    [7]="8007" \
                    [8]="8008" \
                    [9]="8009" \
                    [10]="8010" \
                    [11]="8011" \
                    [12]="8012" \
                    [13]="8013" \
                    [14]="8014" \
                    [15]="8015" \
                    [16]="8016" \
                    [17]="8017" \
                    [18]="8018" \
                    [19]="8019" \
                    [20]="8020" \
                    [21]="8021" \
                    [22]="8022" \
                    [23]="8023" \
                    [24]="8024" \
                    [25]="8025" \
                    [26]="8026" \
                    [27]="8027" \
                    [28]="8028" \
                    [29]="8029" \
                    [30]="8030" \
                    [31]="8031" \
                    [32]="8032" \
                    )


# Use the task array ID to get the right value for this job
# These are defined with the SBATCH header
ptNum=${labels[${SLURM_ARRAY_TASK_ID}]}







#--------- modules ---------

module load matlab/R2019a
module load spm/spm12



#--------- Run ---------

echo; echo; echo; echo;
echo "========== rdm VOI =========="
echo "dir: '${root_dir}' | pt: ${ptNum} | 'name: ${name}' | fwhm: ${fwhm} | confound: '${confound}' | cvi: '${cvi}'"
echo; echo; echo; echo;

echo 'started at:'
date
echo; echo; echo; echo;

# do level 1
matlab-threaded –nodisplay -nodesktop -r "RDM_eigenvariates('${root_dir}','${spm_dir}', ${ptNum}, '${name}', ${fwhm}, '${confound}', '${cvi}');"

echo; echo; echo; echo;
echo 'finished at:'
date


