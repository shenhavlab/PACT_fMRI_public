function RDM_wholeBrain_2_2_contrast(root_dir, spm_dir, varargin)
%% first level analysis - parametric target & distractor
% Harrison Ritz 2021

%
% root_dir = '/Volumes/hritz/data/mri-data/RDM'
% ptNum = 9004

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
pt_dir   = sprintf('%s/spm-data/*/level 1/%s', root_dir, analysis);


%% === get conditions

wt_name = [];
wt_val  = [];
wl = dir(fullfile(pt_dir, 'contrasts.mat'));
load(fullfile(wl(1).folder, wl(1).name));

cons = find(ismember(wt_kind, 'T'));




%% === RUN JOB

spm_jobman('initcfg')


matlabbatch = [];

for cc = 1:length(cons)
        
    
    matlabbatch = [];
    
    matlabbatch{1}.spm.stats.con.spmmat = {sprintf('%s/spm-data/level 2/%s/%s/SPM.mat', root_dir, analysis, wt_name{cons(cc)})};
    
    
    % F
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.name      = 'F';
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights   = 1;
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep   = 'none';
    
    % T +
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name      = 'T+';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights   = 1;
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep   = 'none';
    
    % T -
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name      = 'T-';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights   = -1;
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep   = 'none';
    
    matlabbatch{1}.spm.stats.con.delete = 1;
    
    
    % run
    tic
    spm_jobman('run', matlabbatch);
    toc
    
    
end















end






