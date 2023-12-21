# PACT_fMRI_public
 fMRI and behavioral analysis scripts for Ritz & Shenhav (in press) 'Orthogonal neural encoding of targets and distractors supports multivariate cognitive control'


## 0_behavior

behavior & data files
- data files have all event information 
- analyze_scanner_behav.m: analyze behavior


## 1_xnat2bids

Brown-specific files for getting data on sever + BIDS formatting


## 2_fmriprep

fmriprep preprocessing


## 3_preproc

moving files around & smoothing


## 4_wholeBrain

GLM analysis

- sh_wholeBrain_1_level1.sh: run first-level GLM
- sh_wholeBrain_2_level2.sh: run second-level GLM
- sh_wholeBrain_3_TFCE.sh: TFCE correction
- sh_wholeBrain_4_RSA.sh: just run first-level models for EGA 

- RDM_wholeBrain_1_1_estimateCat.m: core function for first-level analyses


## 5_ROI

Rostral-caudal analysis within dACC ROI


## 6_EGA

EGA analysis

- RDM_parcel_1_1_fit.m: calculate similarity metrics for different models
- RDM_parcel_2_1_analyze.m: analyze similarity metrics
- RDM_parcel_2_3_mediation: run FC analyses


## 7_FC

- extract eigenvariate for FC analyses


## simulations

EGA simulations



## utils

lots of helpful packages

- pcm_toolbox
	- pcm_estGCrossval.m: core similarity metric function

- rsatoolbox
	- basis of parcel-wise searchlight function
	- noise normalization