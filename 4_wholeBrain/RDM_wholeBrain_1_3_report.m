function RDM_targdist_lvl1_report(root_dir, spm_dir, ptNum, varargin)
%% first level contrast - parametric target & distractor
% Harrison Ritz 2021

%
% root_dir = '/Volumes/hritz/data/mri-data/RDM'
% ptNum = 9004



%% === set variables
addpath(genpath(spm_dir)); % add spm12 folder with bug-squashed rwls
spm('defaults', 'fmri')

% set default RAM & use RAM for analysis
spm_get_defaults('maxmem', 128 * 2^30)
spm_get_defaults('resmem', true)
spm_get_defaults('cmdline', true)



%% load variables
if length(varargin) >= 1 && ~isempty(varargin{1})
    name = varargin{1};
else
    name        = 'targdist';
end


if length(varargin) >= 2 && ~isempty(varargin{2})
    fwhm = varargin{2};
else
    fwhm        = 6;
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
participant = ptNum

pt_dir      = sprintf('%s/spm-data/sub-%d', root_dir, ptNum);
save_dir    = sprintf('%s/level 1/%s', pt_dir, analysis);


spmFile = fullfile(save_dir, 'SPM.mat');




%% === run batch
spm_jobman('initcfg')

matlabbatch{1}.spm.stats.results.spmmat = {spmFile};
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{1}.spm.stats.results.conspec.extent = 5;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.export{1}.png = true;



% run job

tic
spm_jobman('run', matlabbatch);
toc











end






