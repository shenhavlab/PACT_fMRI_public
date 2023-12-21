README


Preprocess:

	1) preproc_1_import.sh
		- transfer preprocessed data from BIDS folder to spm-data
		- manually: make sure results files are in RDM/behavior folder

	2) preproc_2_smooth.sh
		- smooth preprocessed data

Whole brain:


	1) wholeBrain_1_level1.sh
		- do first-level analysis
		- estimate --> contrast --> report

	2) wholeBrain_2_level2.sh
		- do second-level analysis
		- estimate --> contrast