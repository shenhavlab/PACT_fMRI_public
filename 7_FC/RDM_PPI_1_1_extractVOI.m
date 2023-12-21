function RDM_PPI_1_1_extractVOI(root_dir, spm_dir, ptNum, varargin)
%% ROI analysis
% Harrison Ritz 2022



voi_masks = {'IPS.nii'}
voi_name = {'IPS'}


%% === set variables
addpath(genpath(spm_dir)); % add spm12 folder with bug-squashed rwls
spm('defaults', 'fmri')

% set default RAM & use RAM for analysis
spm_get_defaults('maxmem', 128 * 2^30)
spm_get_defaults('resmem',true)
spm_get_defaults('cmdline',true)



%% load variables
if length(varargin) >= 1 && ~isempty(varargin{1})
    name = varargin{1};
else
    name        = 'feature';
end

if length(varargin) >= 2 && ~isempty(varargin{2})
    fwhm = varargin{2};
else
    fwhm        = 8;
end

if length(varargin) >= 3 && ~isempty(varargin{3})
    confound    = varargin{3};
else
    confound    = 'movt';
end

if length(varargin) >= 4 && ~isempty(varargin{4})
    cvi         = varargin{4};
else
    cvi         = 'wls';
end




%% folders & names
analysis    = sprintf('%s_s-%dmm_cfd-%s_cvi-%s', name, fwhm, confound, cvi)

pt_dir      = sprintf('%s/spm-data/sub-%d', root_dir, ptNum)
spm_dir     = sprintf('%s/level 1/%s', pt_dir, analysis)
mask_dir    = sprintf('%s/RDM_fmri_scripts/masks', root_dir)
save_dir    = sprintf('%s/spm-data/PPI/data/', root_dir)


% make directories
mkdir(save_dir); % save results

%% get convolved predictor

% s = load(fullfile(spm_dir, 'SPM.mat'));
% 
% dist = s.SPM.xX.X(:, ismember(s.SPM.xX.name, 'Sn(1) trialxdist^1*bf(1)'));
% target = s.SPM.xX.X(:, ismember(s.SPM.xX.name, 'Sn(1) trialxtarget^1*bf(1)'));
% 
% save(fullfile(save_dir, sprintf('task_%d', ptNum)), 'target', 'dist')




%% get VOI ======= re-write for block design

regName = {'targ', 'dist'}


for cc = 1:length(regName)
    for vv = 1:length(voi_masks)


        %% === get VOI
        clear matlabbatch
        matlabbatch{1}.spm.util.voi.spmmat = {fullfile(spm_dir, 'SPM.mat')};
        matlabbatch{1}.spm.util.voi.adjust = nan;
        matlabbatch{1}.spm.util.voi.session = 1;
        matlabbatch{1}.spm.util.voi.name = sprintf('%s_%s_%d', voi_name{vv}, regName{cc}, ptNum);
        matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {fullfile(mask_dir, voi_masks{vv})};
        matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
        matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {fullfile(pt_dir, 'anat', sprintf('sub-%d_ses-01_acq-mprage_space-MNI152NLin6Asym_res-2_desc-brain_mask.nii', ptNum))};
        matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
        matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';


        tic
        spm_jobman('run', matlabbatch);
        toc


        %% === get PPI
        clear matlabbatch
        matlabbatch{1}.spm.stats.ppi.spmmat = {fullfile(spm_dir, 'SPM.mat')};
        matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {fullfile(spm_dir, sprintf('VOI_%s_%s_%d_1.mat', voi_name{vv}, regName{cc}, ptNum))};
        if cc ==1
            matlabbatch{1}.spm.stats.ppi.type.ppi.u = [1 4 1]; % target
        elseif cc == 2
            matlabbatch{1}.spm.stats.ppi.type.ppi.u = [1 5 1]; % distractor
        end
        matlabbatch{1}.spm.stats.ppi.name = sprintf('%s_%s_%d', voi_name{vv}, regName{cc}, ptNum);
        matlabbatch{1}.spm.stats.ppi.disp = 0;


        tic
        spm_jobman('run', matlabbatch);
        toc


        %% move file
        d = dir(sprintf('%s/*%s_%s_%d*', spm_dir, voi_name{vv}, regName{cc}, ptNum))
        copyfile(fullfile(d(1).folder, d(1).name), fullfile(save_dir, d(1).name), 'f');


    end

end


