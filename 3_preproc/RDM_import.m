function RDM_import(root_dir,deriv_dir, ptNum)
% % === gunzip & import images

% root_dir = '/Volumes/hritz/data/mri-data/RDM'
% ptNum = 9004


data_dir = 'spm-data';
cd(root_dir);

mkdir(data_dir);

%% import data

bids_pt_dir = sprintf('%s%s/sub-%d/', root_dir, deriv_dir, ptNum);

%  ====== make directories
pt_rootDir = fullfile(root_dir, data_dir, ['sub-', num2str(ptNum)]);
pt_funcDir = fullfile(pt_rootDir, 'func');
pt_anatDir = fullfile(pt_rootDir, 'anat');
pt_confoundDir = fullfile(pt_rootDir, 'confound');


mkdir(pt_rootDir)
mkdir(pt_funcDir)
mkdir(pt_anatDir)
mkdir(pt_confoundDir)



% %% % ====== FUNC unzip


% load preproc
ft = sprintf('%s/*/func/%s', bids_pt_dir, 'sub*space-MNI152NLin6Asym*desc-preproc_bold.nii.gz')
fl = dir(ft)

fn = [];
for ff = 1:length(fl)
    fn{ff,1} = fullfile(fl(ff).folder, fl(ff).name);
end

disp(fn)



%% % ====== ANAT unzip
ft = sprintf('%s/*/anat/%s', bids_pt_dir, 'sub*space-MNI152NLin6Asym*desc-preproc_T1w.nii.gz')
fl = dir(ft)

fn = [];
for ff = 1:length(fl)
    fn{ff,1} = fullfile(fl(ff).folder, fl(ff).name);
end

disp(fn)

tic
gunzip(fn, pt_anatDir)
toc



%% % ====== BRAIN MASKS unzip
ft = sprintf('%s/*/anat/%s', bids_pt_dir, 'sub*space-MNI152NLin6Asym*desc-brain_mask.nii.gz')
fl = dir(ft);

fn = [];
for ff = 1:length(fl)
    fn{ff,1} = fullfile(fl(ff).folder, fl(ff).name);
end

disp(fn)

tic
gunzip(fn, pt_anatDir)
toc



ft = sprintf('%s/*/anat/%s', bids_pt_dir, 'sub*space-MNI152NLin6Asym*label-GM_probseg.nii.gz')
fl = dir(ft);

fn = [];
for ff = 1:length(fl)
    fn{ff,1} = fullfile(fl(ff).folder, fl(ff).name);
end

disp(fn)

tic
gunzip(fn, pt_anatDir)
toc






%% ====== confound regressors
ft = sprintf('%s/*/func/%s', bids_pt_dir, 'sub*desc-confounds_timeseries.tsv')
fl = dir(ft)

for ff = 1:length(fl)
    copyfile(fullfile(fl(ff).folder, fl(ff).name), pt_confoundDir)
end






end





