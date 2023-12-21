function RDM_wholeBrain_3_1_TFCE(root_dir, spm_dir, varargin)
%% first level analysis - parametric target & distractor
% Harrison Ritz 2021

%
% root_dir = '/Volumes/hritz/data/mri-data/RDM'
% ptNum = 9004

nPerm = 10000

% maskName = 'Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii'
% tfceDir = sprintf('WB__tfce%d', nPerm)

maskName = 'kong22_dACC.nii'
tfceDir = sprintf('dACC__tfce%d', nPerm)



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
pt_dir      = sprintf('%s/spm-data/*/level 1/%s', root_dir, analysis)


%% === get conditions
wt_name = [];

wl = dir(fullfile(pt_dir, 'contrasts.mat'));
load(fullfile(wl(1).folder, wl(1).name));

cons = find(ismember(wt_kind, 'T'));

 %% === get mask
 mask = fullfile(root_dir, 'RDM_fmri_scripts', 'masks', maskName)



%% === estimate second level for each condition
spm_jobman('initcfg')

for cc = 1:length(cons)
    
    % define save path
    dirPALM = sprintf('%s/RDM_fmri_scripts/util/palm', root_dir)
    save_dir = sprintf('%s/spm-data/level 2/%s/%s/%s', root_dir, analysis, wt_name{cons(cc)}, tfceDir);
    dirIN = fullfile(save_dir, 'con4D.nii')
    dirOUT = fullfile(save_dir, 'T')
    
    % make save path
    if exist(save_dir, 'dir')
        rmdir(save_dir, 's'); % remove existing results folder
    end
    mkdir(save_dir);
    
    


    % get files
    cl = dir(fullfile(pt_dir, sprintf('con_%.4d.nii', cons(cc))));
    cf = [];
    for ci = 1:length(cl)

        % read
        cf = fullfile(cl(ci).folder, cl(ci).name);
        VolIn = spm_vol(cf);

        % prep
        Vo      = struct(...
            'fname',    dirIN,...
            'dim',      VolIn(1).dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      VolIn(1).mat,...
            'n',        [ci 1],...
            'descrip',  'con');

        % write
        spm_write_vol(Vo, spm_read_vols(VolIn));

    end

    
    
    % batch -- convert to 4d
    %     matlabbatch = [];
    %
    %     matlabbatch{1}.spm.util.cat.vols = cf;
    %     matlabbatch{1}.spm.util.cat.name = dirIN;
    %     matlabbatch{1}.spm.util.cat.dtype = 0;
    %     matlabbatch{1}.spm.util.cat.RT = NaN;
    %
    %
    %     % === RUN BATCH
    %     tic
    %     spm_jobman('run', matlabbatch);
    %     toc
    
    
    % === TFCE =============
    system(sprintf('sbatch --export=dirPALM="%s",mask="%s",dirIN="%s",dirOUT="%s",nPerm=%d runTFCE.sh', dirPALM, mask, dirIN, dirOUT, nPerm))
    
    
end

end