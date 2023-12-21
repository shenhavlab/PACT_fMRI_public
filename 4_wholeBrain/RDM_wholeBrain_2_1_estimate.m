function RDM_wholeBrain_2_1_estimate(root_dir, spm_dir, varargin)
%% first level analysis - parametric target & distractor
% Harrison Ritz 2021


% save out mask
maskName = 'kong22_dACC.nii'


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
pt_dir   = sprintf('%s/spm-data/*/level 1/%s', root_dir, analysis);


%% === get conditions
wt_name = [];

wl = dir(fullfile(pt_dir, 'contrasts.mat'));
load(fullfile(wl(1).folder, wl(1).name));

cons = find(ismember(wt_kind, 'T'));


%% get mask

addpath(genpath(sprintf('%s/RDM_fmri_scripts/util/NIfTI_Toolbox', root_dir)))
mask = load_nii(fullfile(sprintf('%s/RDM_fmri_scripts/masks', root_dir), maskName));

[mX, mY, mZ] = ind2sub(size(mask.img), find(mask.img));
mXYZ = [mX, mY, mZ]';




%% === estimate second level for each condition
spm_jobman('initcfg')


% new analysis dir
analysis_dir = sprintf('%s/spm-data/level 2/%s/%s', root_dir, analysis);
if exist(analysis_dir, 'dir')
    rmdir(analysis_dir, 's');
end
mkdir(analysis_dir);

clear data
for cc = 1:length(cons)
    
    % make save path
    save_dir = fullfile(analysis_dir, wt_name{cons(cc)}); % sprintf('%s/spm-data/level 2/%s/%s', root_dir, analysis, wt_name{cons(cc)});    
    
    
    % get files
    cl = dir(fullfile(pt_dir, sprintf('con_%.4d.nii', cons(cc))));
    cf = [];
    for ci = 1:length(cl)
        cf{ci,1} = fullfile(cl(ci).folder, cl(ci).name);
    end
        
    
    % batch
    matlabbatch = [];
    
    % specify
    matlabbatch{1}.spm.stats.factorial_design.dir = {save_dir};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cf;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % estimate
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    
    % === RUN
    tic
    spm_jobman('run', matlabbatch);
    toc


    % === save masked data
    data{cc} = spm_get_data(cf, mXYZ);
    
end

save(fullfile(analysis_dir, maskName(1:end-4)), 'data', 'mXYZ')



end