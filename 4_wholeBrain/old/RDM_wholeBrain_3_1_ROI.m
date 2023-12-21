function RDM_wholeBrain_3_1_ROI(root_dir, spm_dir, varargin)
%% ROI analysis
% Harrison Ritz 2022

%
% root_dir = '/Volumes/hritz/data/mri-data/RDM'
% ptNum = 9004

conList = {'dist', 'rt', 'acc', 'distRT', 'distAcc', 'rtAcc', 'distRtAcc'}

mask = 'Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii'

maskList = [285, 286, 287, 295,...
    255,...
    281, 282, 283, 293, 294, ...
    85, 86, 90, ...
    73,...
    80, 81,...
    274,...
    302, 303,...
    56,...
    79,...
    99,100, 106]



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
    name        = 'targdist';
end

if length(varargin) >= 2 && ~isempty(varargin{2})
    fwhm = varargin{2};
else
    fwhm        = 8;
end

if length(varargin) >= 3 && ~isempty(varargin{3})
    confound    = varargin{3};
else
    confound    = 'full';
end

if length(varargin) >= 4 && ~isempty(varargin{4})
    cvi         = varargin{4};
else
    cvi         = 'wls';
end



%% folders & names
analysis    = sprintf('%s_s-%dmm_cfd-%s_cvi-%s', name, fwhm, confound, cvi)
l2_dir      = sprintf('%s/spm-data/level 2/%s', root_dir, analysis)
save_dir    = sprintf('%s/RDM_fmri_scripts/4_wholeBrain/ROI/%s', root_dir, analysis)
mask_dir    = sprintf('%s/RDM_fmri_scripts/masks', root_dir)


% make save path
if exist(save_dir, 'dir')
    rmdir(save_dir, 's'); % remove existing results folder
end
mkdir(save_dir);

%% extract betas
clear D P T

% make brain template
V = spm_vol(fullfile(mask_dir, mask));
Y = spm_read_vols(V,1);
blank_brain = zeros(size(Y));


for cc = 1:length(conList) % for each contrast


    % make plotting ROI
    [V_T, V_P] = deal(V);

    [brain_T, brain_P] = deal(blank_brain);


    [V_T.fname, V_T.private.dat.fname] = deal(fullfile(save_dir, [conList{cc}, '_T.nii']));
    V_T.descrip = sprintf('ROI analysis: %s tval', conList{cc});

    [V_P.fname, V_P.private.dat.fname] = deal(fullfile(save_dir, [conList{cc}, '_logP.nii']));
    V_P.descrip = sprintf('ROI analysis: %s pval', conList{cc});


    % load SPM
    s = load(fullfile(l2_dir, conList{cc}, 'SPM.mat'));


    for rr = 1:length(maskList) % for each ROI


        % load ROI
        indx = find(Y ==  maskList(rr));
        [x,y,z] = ind2sub(size(Y),indx);

        XYZ = [x y z]';


        % data
        bold = spm_get_data(s.SPM.xY.P, XYZ);
        D{rr,cc} = bold;

        % stats
        [~, pval, ~, stats] = ttest(nanmean(bold,2)./nanstd(bold,[],2));
        P(rr,cc) = pval;
        T(rr,cc) = stats.tstat;

        % plot ROI
        brain_T(indx) = stats.tstat;
        brain_P(indx) = abs(log10(pval));


    end

    % write ROI map
    spm_write_vol(V_T, brain_T)
    spm_write_vol(V_P, brain_P)

end


%% targ-dist correlation
clear targdist_R
if sum(ismember(conList, {'dist', 'target'}))==2


    % make plotting ROI
    [V_R, V_RlogP] = deal(V);
    [brain_R, brain_RlogP] = deal(blank_brain);

    [V_R.fname, V_R.private.dat.fname] = deal(fullfile(save_dir, 'targdist_R.nii'));
    V_R.descrip = sprintf('ROI analysis: %s rval', conList{cc})

    [V_RlogP.fname, V_RlogP.private.dat.fname] = deal(fullfile(save_dir, 'targdist_RlogP.nii'));
    V_RlogP.descrip = sprintf('ROI analysis: %s r pval', conList{cc})


    for rr = 1:length(maskList) % for each ROI

        % r val
        targdist_R(rr,cc) = mean(diag(corr(D{rr, ismember(conList, 'dist')}', D{rr,ismember(conList, 'target')}', 'rows', 'pairwise')));
        bf10 = bf_ttest(atanh(diag(corr(D{rr, ismember(conList, 'dist')}', D{rr,ismember(conList, 'target')}', 'rows', 'pairwise'))));

        % write
        brain_R(Y ==  maskList(rr)) = targdist_R(rr,cc);
        brain_RlogP(Y ==  maskList(rr)) = log10(bf10);

    end

    spm_write_vol(V_R, brain_R)
    spm_write_vol(V_RlogP, brain_RlogP)

else

    targdist_R = nan;

end





%% targ-dist correlation
clear targdist2_R
if sum(ismember(conList, {'dist2', 'target'}))==2


    % make plotting ROI
    [V_R, V_RlogP] = deal(V);
    [brain_R, brain_RlogP] = deal(blank_brain);

    [V_R.fname, V_R.private.dat.fname] = deal(fullfile(save_dir, 'targdist2_R.nii'));
    V_R.descrip = sprintf('ROI analysis: %s rval', conList{cc})

    [V_RlogP.fname, V_RlogP.private.dat.fname] = deal(fullfile(save_dir, 'targdist2_RlogP.nii'));
    V_RlogP.descrip = sprintf('ROI analysis: %s r pval', conList{cc})


    for rr = 1:length(maskList) % for each ROI

        % r val
        targdist2_R(rr,cc) = mean(diag(corr(D{rr, ismember(conList, 'dist2')}', D{rr,ismember(conList, 'target')}', 'rows', 'pairwise')));
        bf10 = bf_ttest(atanh(diag(corr(D{rr, ismember(conList, 'dist2')}', D{rr,ismember(conList, 'target')}', 'rows', 'pairwise'))));

        % write
        brain_R(Y ==  maskList(rr)) = targdist2_R(rr,cc);
        brain_RlogP(Y ==  maskList(rr)) = log10(bf10);

    end

    spm_write_vol(V_R, brain_R)
    spm_write_vol(V_RlogP, brain_RlogP)

else

    targdist2_R = nan;

end







%% save stats
save(fullfile(save_dir, 'stats'), 'conList', 'maskList', 'D', 'P', 'T', 'targdist_R', 'targdist2_R', '-v7.3')






