function RDM_wholeBrain_1_2_contrast(root_dir, spm_dir, ptNum, varargin)
%% second level analysis ===============================================
% Harrison Ritz 2021



%% === initialize
addpath(genpath(spm_dir)); % add spm12 folder
spm('defaults', 'fmri')


% set default RAM & use RAM for analysis
spm_get_defaults('maxmem',64 * 2^30)
spm_get_defaults('resmem',true)
spm_get_defaults('cmdline',true)


%% load variables
if length(varargin) >= 1 && ~isempty(varargin{1})
    name = varargin{1};
else
    name        = 'linear';
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

pt_dir      = sprintf('%s/spm-data/sub-%d', root_dir, ptNum);
save_dir    = sprintf('%s/level 1/%s', pt_dir, analysis);
cond_dir    = sprintf('%s/conditions', pt_dir);

spmFile = fullfile(save_dir, 'SPM.mat');



%% === load confounds
load(sprintf('%s/%s_condList.mat', cond_dir, analysis));

%% === set contrasts




[wt_name, wt_val, wt_kind] = deal([]);



switch name

    
  case 'feature' % =========================================================
        %%

        conNames = {...
            'trial', 'targResp', 'distResp', 'targCoh', 'distCoh', 'cong', 'lapse'};

        findCond = @(name) (condList == find(ismember(conNames, name))) ./ sum(condList == find(ismember(conNames, name)));


        nc = 1;

        for cc = 2:(length(conNames)-1)

            wt_name{nc} = conNames{cc};
            wt_val{nc}  = findCond(conNames{cc});
            wt_kind{nc} = 'T';
            nc = nc+1;

        end


  



end











save(fullfile(save_dir, 'contrasts'), 'wt_name', 'wt_val', 'wt_kind');




%% === RUN JOB
spm_jobman('initcfg')


matlabbatch = [];

matlabbatch{1}.spm.stats.con.spmmat = {spmFile};


for ww = 1:length(wt_name)

    switch wt_kind{ww}

        case 'F'

            matlabbatch{1}.spm.stats.con.consess{ww}.fcon.name      = wt_name{ww};
            matlabbatch{1}.spm.stats.con.consess{ww}.fcon.weights   = double(wt_val{ww});
            matlabbatch{1}.spm.stats.con.consess{ww}.fcon.sessrep   = 'none';

        case 'T'

            matlabbatch{1}.spm.stats.con.consess{ww}.tcon.name      = wt_name{ww};
            matlabbatch{1}.spm.stats.con.consess{ww}.tcon.weights   = double(wt_val{ww});
            matlabbatch{1}.spm.stats.con.consess{ww}.tcon.sessrep   = 'none';

    end

end

matlabbatch{1}.spm.stats.con.delete = 1;






% run
tic
spm_jobman('run', matlabbatch);
toc











end






