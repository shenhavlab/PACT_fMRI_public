function RDM_parcel_1_1_fit(root_dir, spm_dir, ptNum, varargin)
%% first level analysis ===================================================
% Harrison Ritz 2021


%% load arguments
if length(varargin) >= 1 && ~isempty(varargin{1})
    name = varargin{1};
else
    name        = 'rdm';
end


if length(varargin) >= 2 && ~isempty(varargin{2})
    fwhm = varargin{2};
else
    fwhm        = 0;
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



%% high-level parameters
name

switch name


    case 'featureBlk'
        %%
        analysisName    = 'feature_parcel'
        condLabel       = {'respTarg', 'respDist', 'cohTarg', 'cohDist', 'cong'}
        condVec         = [1:length(condLabel)]'

    case 'FC_PFCl'
        %%
        analysisName    = 'FC_PFCl_parcel'
        condLabel       = {...
            'respTarg', 'respDist', 'cohTarg',  'cohDist', 'cong', 'rt', 'acc', 'PFCl' ...
            };
        condVec         = [1:length(condLabel)-1]'


    case 'FC_IPS-PFCl'
        %%
        analysisName    = 'FC_IPS-PFCl_parcel'
        condLabel       = {...
            'respTarg', 'respDist', 'cohTarg',  'cohDist', 'cong', 'rt', 'acc', 'IPS', 'PFCl' ...
            };
        condVec         = [1:length(condLabel)-2]'


    case 'FC_IPS'
        %%
        analysisName    = 'FC_IPS_parcel'
        condLabel       = {...
            'respTarg', 'respDist', 'cohTarg',  'cohDist', 'cong', 'rt', 'acc', 'IPS' ...
            };
        condVec         = [1:length(condLabel)-1]'


    case 'perfCong'
        %%
        analysisName    = 'perfCong'
        condLabel       = {...
            'respTarg', 'respDist', 'cohTarg',  'cohDist', 'cong', 'rt', 'acc', ...
            };
        condVec         = [1:length(condLabel)]'


    case 'perfFullBlk-comb'
        %%
        analysisName    = 'perfFullBlkTask_parcel'
        condLabel       = {...
            'motRespTarg', 'motRespDist', 'motCohTarg', 'motCohDist', 'motCong', ...
            'colRespTarg', 'colRespDist', 'colCohTarg', 'colCohDist', 'colCong'};
        condVec         = [1:length(condLabel)]'

    case 'perfCongBlk'
        %%
        analysisName    = 'perf_parcel'
        condLabel       = {'respTarg', 'respDist', 'cohTarg',  'cohDist', 'cong', 'rt', 'acc'}
        condVec         = [1:length(condLabel)]'


    case 'perfBlk-comb'
        %%
        analysisName    = 'perfBlkTask_parcel'
        condLabel       = {...
            'motRespTarg', 'motRespDist', 'motCohTarg', 'motCohDist', 'motCong', ...
            'colRespTarg', 'colRespDist', 'colCohTarg', 'colCohDist', 'colCong'};
        condVec         = [1:length(condLabel)]'


    case 'featureBlk-comb'
        %%
        analysisName    = 'featureBlkTask_parcel'
        condLabel       = {'motRespTarg', 'motRespDist', 'motCohTarg', 'motCohDist', 'motCong', 'colRespTarg', 'colRespDist', 'colCohTarg',  'colCohDist', 'colCong'}
        condVec         = [1:length(condLabel)]'




    case 'margRespBlk'
        %%
        analysisName    = 'margResp_parcel'
        condLabel       = {'targ4R', 'targ3R', 'targ2R',  'targ1R', 'targ1L', 'targ2L', 'targ3L',  'targ4L', 'dist2R', 'dist1R', 'dist1L',  'dist2L'}
        condVec         = [1:length(condLabel)]'


    case 'perfBlk'
        %%
        analysisName    = 'perf_parcel'
        condLabel       = {'respTarg', 'respDist', 'cohTarg',  'cohDist', 'rt', 'acc'}
        condVec         = [1:length(condLabel)]'



end



%% === initialize
addpath(genpath(spm_dir)); % add spm12 folder
spm('defaults', 'fmri')

addpath(genpath(fullfile(root_dir,'RDM_fmri_scripts', 'util')));
import rsa.*


delete('./temp*.mat') % clear temps

% set default RAM & use RAM for analysis
spm_get_defaults('maxmem',64 * 2^30)
spm_get_defaults('resmem',true)
spm_get_defaults('cmdline',true)



%% folders & names
analysis    = sprintf('%s_s-%dmm_cfd-%s_cvi-%s', name, fwhm, confound, cvi)

pt_dir      = sprintf('spm-data/sub-%d', ptNum)
l1_dir      = sprintf('%s/level 1/%s', pt_dir, analysis)
loc_dir     =  sprintf('%s/level 1/%s', pt_dir, sprintf('%s_s-%dmm_cfd-%s_cvi-%s', 'localizer', fwhm, confound, cvi))

mask_dir    = sprintf('%s/RDM_fmri_scripts/masks', root_dir)




% make save_dir
mkdir(save_dir);



%% === get parcellation
parcel.Vo   = spm_vol(fullfile(mask_dir, 'Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii'));
parcel.info = readtable(fullfile(mask_dir, 'Schaefer2018_400Parcels_Kong2022_17Networks_order.lut.txt'));
parcel

%% === conditions & parcels

conds = 0; % skip trial main effect
parts = 1; % skip trial main effect

% -=---------------- UPDATE!
regName = regexp(name, {'Half'});
if any(regName{1})
    fprintf('\n\nHALF\n\n')
    nBlk = 2;
else
    nBlk = (6 - ismember(ptNum, [8010, 8019]));
end


for ii = 1:nBlk

    conds = [conds; condVec];
    parts = [parts; ii*ones(size(condVec))];

end


switch name


    case 'FC_M1Left'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FC_VISA-IPS'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FC_VISA'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FC_COL-IPS'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FC_COL'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FC_PFCl'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FC_PFCl_Sal'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FC_PFCl_Ctrl'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FC_MOT-IPS'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FC_MOT-PFCl'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FC_MOT'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FC_IPS-PFCl'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FC_IPS-PFCl_Sal'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];


    case 'FC_IPS-PFCl_Ctrl'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];


    case 'FC_IPS-dACC'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FC_IPS'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];


    case 'FC_SPL'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];


    case 'FC_dACC-sal-PFCl'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'shift_dACC-sal'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];


    case 'PPI_PFCl'
        conds = [conds; 0;...
            repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1]); repmat(condVec(end)+3, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'; [1:nBlk]'];


    case 'FC_dACC-sal'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FC_dACC-ctrl'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'PPI_dACC-subnet'
        conds = [conds; 0;...
            repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1]); repmat(condVec(end)+3, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'; [1:nBlk]'];

    case 'featurePPI'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FC_dACC' % add FC
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FC_dACC-HG' % add FC
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FC_dACC-subnet' % add FC
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FC_dACC-pc' % add FC
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FC_dACC-PFCl' % add FC
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FCperf_dACC-subnet' % add FC
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FCperf2_dACC-subnet'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FCperfD_dACC-subnet'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FCperf_dACC-PFCl'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'];

    case 'FCperf_dACC'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FCperf_PFCl'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'];

    case 'FCperf_dACC-subnet-PFCl'
        conds = [conds; 0; repmat(condVec(end)+1, [nBlk, 1]); repmat(condVec(end)+2, [nBlk, 1]); repmat(condVec(end)+3, [nBlk, 1])];
        parts = [parts; 1; [1:nBlk]'; [1:nBlk]'; [1:nBlk]'];

end


Opt.analysisName    = analysisName;     %  Name of the analysis
Opt.rootPath        = root_dir;
Opt.spmDir          = l1_dir;           %  Directory that contains the first-level analysis
Opt.conditionVec    = conds;            %  Vector, indicating which of the betas belong to which condition
Opt.conditionLabels = {};               %  Condition labels - used if conditionVec is not given
Opt.partition       = parts;            %  Indicator of partition
Opt.imageDataDir    = [];               %  Directory where the pre-processed raw data resides (if different from what is specified in the SPM-structure)
Opt.saveSigma       = 'none';           %  Determine whether to save the variance / covariance of the beta estimates as well as the distances
Opt.saveFn          = fullfile(save_dir, sprintf('sub%d_%s', ptNum, analysisName));
Opt.nBlk            = nBlk;
Opt.condLabel       = condLabel;
Opt.loc_dir         = loc_dir;
Opt


%% === run RSA

disp('starting fit')
wholeTime = tic;


runSearchlightParcel(parcel, Opt)


toc(wholeTime)
disp('finished fit, good job!')



