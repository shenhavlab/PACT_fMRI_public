#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --mem=128GB
#SBATCH -N 1
#SBATCH -n 18
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=harrison_ritz@brown.edu
#SBATCH --account=carney-ashenhav-condo  
#SBATCH -J RDM-fmriprep
#SBATCH --output logs/fmriprep-log_%A-%a.txt
#SBATCH --error logs/fmriprep-error_%A-%a.txt
#SBATCH --array=4-32

#--------- CONFIGURE THESE VARIABLES ---------

bids_root_dir=/gpfs/data/ashenhav/mri-data/RDM2     # based on oscar path
investigator=ashenhav                               # investagor   
study_label=1222                                    # study label 
fmriprep_version=20.2.6                             # check /gpfs/data/bnc/simgs/poldracklab for the latest version
nthreads=18                                         # heuristic: number of cores
nDummy=2                                            # number of dummy scans


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
SUBJ_LABEL=${labels[${SLURM_ARRAY_TASK_ID}]}

#--------- Run FMRIPREP --------- 

# runs singularity for each subject
singularity run --cleanenv                                              	        \
    --bind ${bids_root_dir}/${investigator}/study-${study_label}:/data              \
    --bind /gpfs/scratch/hritz:/scratch                                             \
    --bind /gpfs/data/bnc/licenses:/licenses                                        \
    /gpfs/data/bnc/simgs/nipreps/fmriprep-${fmriprep_version}.sif                   \
    /data/bids                                                                      \
    /data/bids/derivatives/fmriprep-${fmriprep_version}                             \
    participant 															        \
    --participant-label ${SUBJ_LABEL} 								                \
    --fs-license-file /licenses/freesurfer-license.txt 						        \
    -w /scratch/fmriprep/${SUBJ_LABEL} 												\
    --output-spaces MNI152NLin6Asym:res-2                                           \
    --nthreads ${nthreads} 															\
    --write-graph 															        \
    --dummy-scans ${nDummy}                                                         \
    --stop-on-first-crash                                                           \
    # --use-aroma                                                                     \
    # --error-on-aroma-warnings                                                       \
    # --fs-no-reconall                                                                \





