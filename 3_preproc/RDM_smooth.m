function RDM_smooth(root_dir, funcRegex, ptNum, fwhm)
%% smooth functional images
% Harrison Ritz 2020

% root_dir = '/Volumes/hritz/data/mri-data/RDM'
% ptNum = 9004

%% === parameters 
spm('defaults', 'fmri')

% set default RAM & use RAM for analysis
spm_get_defaults('maxmem',64 * 2^30)
spm_get_defaults('resmem',true)
spm_get_defaults('cmdline',true)


% parameters
smooth_prefix   = sprintf('s-%dmm_', fwhm)
func_dir        = sprintf('%s/spm-data/sub-%d/func/%s', root_dir, ptNum)


% clear existing smoothed data
delete(fullfile(func_dir, sprintf('%s*', smooth_prefix)));



% === get files

fl = dir(fullfile(func_dir, funcRegex))

fn = [];
for ff = 1:length(fl)
    fn{ff,1} = fullfile(fl(ff).folder, fl(ff).name);
end


%% === BATCH ====
%  ==============

% === setup batch
spm_jobman('initcfg')


% get frames
matlabbatch{1}.spm.util.exp_frames.files = fn;
matlabbatch{1}.spm.util.exp_frames.frames = Inf;

% smooth
matlabbatch{2}.spm.spatial.smooth.data(1) = ...
    cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 1;
matlabbatch{2}.spm.spatial.smooth.prefix = smooth_prefix;



% === run
spm_jobman('run', matlabbatch);



end