#!/bin/bash
#SBATCH -t 04:00:00
#SBATCH --mem=128GB
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hritz@brown.edu
#SBATCH --account=carney-ashenhav-condo  
#SBATCH -J import_pt
#SBATCH --output logs/import_%A-%a.txt
#SBATCH --array=4-32




#--------- Variables ---------
root_dir="/users/hritz/data/mri-data/RDM2"
deriv_dir="/ashenhav/study-1222/bids/derivatives/fmriprep-20.2.6/fmriprep"



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


#--------- Modules ---------

module load matlab/R2019a


#--------- Run ---------

echo; echo; echo; echo;
echo "========== transferring data =========="
echo; echo; echo; echo;

echo 'started at:'
date
echo; echo; echo; echo;



# do smoothing
matlab-threaded –nodisplay -r "RDM_import('${root_dir}', '${deriv_dir}', ${ptNum}); exit"



echo; echo; echo; echo;
echo 'finished at:'
date