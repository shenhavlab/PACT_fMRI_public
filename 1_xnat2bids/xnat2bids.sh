#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH --mem=32GB
#SBATCH -m arbitrary
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --account=carney-ashenhav-condo  
#SBATCH -J RDM-xnat2bids
#SBATCH --output logs/xnat2bids_%A-%a.txt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=harrison_ritz@brown.edu
#SBATCH --array=32



# Reference:
# https://github.com/brown-bnc/sanes_sadlum/blob/main/preprocessing/xnat2bids/run_xnat2bids.sh


#----------- Study variables -----
investigator_name="ashenhav"
study_label=1222

run_xnat2bids=true
run_bidspostprocess=true


#--------- Variables ---------
# This line makes our bash script complaint if we have undefined variables
set -u


# Read variables in the .env file in current directory
# This will read:
# XNAT_USER, XNAT_PASSWORD
# set -a
# [ -f .bashrc ] && . .bashrc
# set +a

# Create temporary token for XNAT, so your password is not saved anywhere!
# To generate, run the following command: bash /gpfs/data/bnc/scripts/xnat-token


XNAT_USER=hritz 					# xnat username
XNAT_PASSWORD=${xnat_password} 		# saved xnat password

# uncomment to print environment variables with XNAT in the name
# printenv | grep XNAT

#--------- xnat-tools ---------
# version of xnat2bids being used
version=v1.0.5        #v1.0.3
# Path to Singularity Image for xnat-tools (maintained by bnc)
simg=/gpfs/data/bnc/simgs/brownbnc/xnat-tools-${version}.sif

#--------- directories ---------
# Your working directory in Oscar, usually /gpfs/data/<your PI's group>.
# We pass (bind) this path to singularity so that it can access/see it
data_dir=/gpfs/data/ashenhav

# Output directory
# It has to be under the data_dir, otherwise it won't be seen by singularity
# bids_root_dir=${data_dir}/shared/bids-export/${USER}
bids_root_dir=${data_dir}/mri-data/RDM2   
mkdir -p -m 775 ${bids_root_dir} || echo "Output directory already exists"

# Bidsmap file for your study
# It has to be under the data_dir, otherwise it won't be seen by singularity
#bidsmap_file=${data_dir}/shared/xnat-tools-examples/${USER}/bidsmaps/sanes_sadlum.json
bidsmap_file=${bids_root_dir}/RDM_fmri_scripts/bidsmaps/bidsmaps_RDM.json	

#----------- Dictionaries for subject specific variables -----
# Check for XNAT Ascension Number Here: https://xnat.bnc.brown.edu/data/experiments/
# Dictionary of sessions per subject
declare -A sessions=(\
					[4]="XNAT16_E00005" \
					[5]="XNAT16_E00004" \
					[6]="XNAT17_E00005" \
					[7]="XNAT19_E00002" \
					[8]="XNAT19_E00001" \
					[9]="XNAT27_E00034" \
					[10]="XNAT27_E00033" \
					[11]="XNAT27_E00038" \
					[12]="XNAT27_E00039" \
					[13]="XNAT27_E00040" \
					[14]="XNAT27_E00041" \
					[15]="XNAT27_E00042" \
					[16]="XNAT27_E00045" \
					[17]="XNAT27_E00046" \
					[18]="XNAT27_E00047" \
					[19]="XNAT27_E00049" \
					[20]="XNAT27_E00048" \
					[21]="XNAT29_E00022" \
					[22]="XNAT29_E00024" \
					[23]="XNAT29_E00026" \
					[24]="XNAT29_E00034" \
					[25]="XNAT29_E00033" \
					[26]="XNAT29_E00036" \
					[27]="XNAT29_E00035" \
					[28]="XNAT29_E00041" \
					[29]="XNAT29_E00040" \
					[30]="XNAT29_E00038" \
					[31]="XNAT29_E00039" \
					[32]="XNAT29_E00046" \
					)



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



# Dictionary of series to skip per subject
declare -A skip_map=(\
					[4]="-s 7" \
					[5]="" \
					[6]="-s 20" \
					[7]="" \
					[8]="" \
					[9]="" \
					[10]="-s 16" \
					[11]="" \
					[12]="" \
					[13]="" \
					[14]="" \
					[15]="" \
					[16]="" \
					[17]="" \
					[18]="" \
					[19]="-s 16" \
					[20]="" \
					[21]="" \
					[22]="" \
					[23]="" \
					[24]="" \
					[25]="" \
					[26]="" \
					[27]="" \
					[28]="" \
					[29]="" \
					[30]="" \
					[31]="" \
					[32]="" \
					)



# Use the task array ID to get the right value for this job
# These are defined with the SBATCH header
echo ${SLURM_ARRAY_TASK_ID}
XNAT_SESSION=${sessions[${SLURM_ARRAY_TASK_ID}]}
SUBJ_LABEL=${labels[${SLURM_ARRAY_TASK_ID}]}
SKIP_STRING=${skip_map[${SLURM_ARRAY_TASK_ID}]}

echo "Processing session:"
echo ${XNAT_SESSION}
echo "Subject label:"
echo ${SUBJ_LABEL}
echo "Series to skip:"
echo ${SKIP_STRING}


#--------- Run xnat2bids ---------
# runs singularity command to extract DICOMs from xnat and export to BIDS
# this command tells singularity to launch out xnat-tools-${version}.sif image
# and execute the xnat2bids command with the given inputs.
# The `-B ${data_dir}` makes that directory available to the singularity container
# The file system inside your container is not the same as in Oscar, unless you bind the paths
# The -i passes a sequence to download, without any -i all sequences will be processed
if "${run_xnat2bids}"; then

 	singularity exec --no-home -B ${data_dir} ${simg} 	\
	 	xnat2bids ${XNAT_SESSION} ${bids_root_dir} 		\
    	-u ${XNAT_USER} 	\
    	-p ${XNAT_PASSWORD} \
    	-f ${bidsmap_file} 	\
		-v --overwrite		\
		${SKIP_STRING} 	
	   	
fi

# -h https://xnat.bnc.brown.edu \

#--------- Run bids-postprocess ---------
# runs singularity command to add "IntendedFor" argument to fieldmap jsons
# NOTE: You will have to verify that the participant_labels are correct above
if "${run_bidspostprocess}"; then

	bids_sub_dir=${bids_root_dir}/${investigator_name}/study-${study_label}/bids
	
	echo "EXPERIMENT DIRECTORY: ${bids_sub_dir}"

	singularity exec --no-home -B ${data_dir} ${simg} 		\
		bids-postprocess ${XNAT_SESSION} ${bids_sub_dir} 	\
		-v -v											\
		--includesubj ${SUBJ_LABEL} 

fi
