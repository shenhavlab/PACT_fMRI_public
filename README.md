fMRI and behavioral analysis scripts for Ritz & Shenhav (in press) 'Orthogonal neural encoding of targets and distractors supports multivariate cognitive control'.

Includes code for Ecoding Geometry Analysis (EGA) and related multivariate functional connectivity.


## 0_behavior

behavior & data files
- data files have all event information 
- `analyze_scanner_behav.m`: analyze behavior


## 1_xnat2bids

Brown-specific files for getting data on sever + BIDS formatting


## 2_fmriprep

fmriprep preprocessing


## 3_preproc

moving files around & smoothing


## 4_wholeBrain

GLM analysis

### batch scripts
- `sh_wholeBrain_1_level1.sh`: run first-level GLM
- `sh_wholeBrain_2_level2.sh`: run second-level GLM
- `sh_wholeBrain_3_TFCE.sh`: TFCE correction
- `sh_wholeBrain_4_RSA.sh`: just run first-level models for EGA 

### functions
- `RDM_wholeBrain_1_1_estimateCat.m`: core function for first-level analyses


## 5_ROI

Rostral-caudal analysis within dACC ROI


## 6_EGA

EGA analysis. First level models need to be fit block-wise for CV similarity.
draws heavily from `pcm_toolbox` & especially `rsatoolbox` (see `utils`)

### batch scripts
- `sh_parcel_1_fit.sh`: run EGA

### functions
- `RDM_parcel_1_1_fit.m`: calculate similarity metrics for different models
- `RDM_parcel_2_1_analyze.m`: analyze similarity metrics
- `RDM_parcel_2_3_mediation`: run FC analyses



## 7_FC

- extract eigenvariate for FC analyses (drawing heavily from SPM's PPI)


## simulations

EGA simulations



## utils

lots of helpful packages

- `\pcm_toolbox\`: 
	- `pcm_estGCrossval.m`: core similarity metric function, adapted with more metrics

- `\rsatoolbox\`
	- basis of parcel-wise searchlight function
	- noise normalization



- `bayesFactor`: https://github.com/klabhub/bayesFactor
- `palm`: https://github.com/andersonwinkler/PALM
- `pcm_toolbox`: https://github.com/jdiedrichsen/pcm_toolbox
- `rsatoolbox`: https://github.com/rsagroup/rsatoolbox_matlab/tree/develop
- `ScientificColourMaps7`: https://www.fabiocrameri.ch/colourmaps/
- `NIfTI_Toolbox`: https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
- `rwls`: https://github.com/jdiedrichsen/rwls
- `surfing`: https://github.com/nno/surfing