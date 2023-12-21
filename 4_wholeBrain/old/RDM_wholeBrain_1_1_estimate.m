function RDM_wholeBrain_1_1_estimate(root_dir, spm_dir, ptNum, varargin)
%% first level analysis - parametric target & distractor
% Harrison Ritz 2021

%
% root_dir = '/Volumes/hritz/data/mri-data/RDM'
% ptNum = 9004


%% === initialize
addpath(genpath(spm_dir)); % add spm12 folder with bug-squashed rwls
spm('defaults', 'fmri')

% set default RAM & use RAM for analysis
spm_get_defaults('maxmem',128 * 2^30)
spm_get_defaults('resmem',true)
spm_get_defaults('cmdline',true)


%% load variables
if length(varargin) >= 1 && ~isempty(varargin{1})
    name = varargin{1};
else
    name        = 'rdm';
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

behav_dir   = sprintf('%s/behavior', root_dir)
pt_dir      = sprintf('%s/spm-data/sub-%d', root_dir, ptNum)
func_dir    = sprintf('%s/func', pt_dir)
save_dir    = sprintf('%s/level 1/%s', pt_dir, analysis)
cond_dir    = sprintf('%s/conditions', pt_dir)
cfd_dir     = sprintf('%s/confound', pt_dir)
sense_dir   = fullfile(root_dir, 'RDM_fmri_scripts', 'behavior')

mask_filt   = 'sub*desc-brain_mask.nii'   % mask filter

regName = regexp(name, {'mot', 'comb'});

if any(regName{1}) % if MOT
    
    movtName = 'motion';
    task_runs   = 1:2:12;
    switch ptNum
        case 8010
            task_runs(task_runs==10) = [];
        case 8019
            task_runs(task_runs==10) = [];
    end
    
    cl = dir(sprintf('%s/sub*task-RDMmotion*.tsv', cfd_dir))
    fl = dir(fullfile(func_dir, sprintf('s-%dmm_*task-RDMmotion*.nii', fwhm)))
    
elseif any(regName{2}) % if COMBINED
    
    movtName = 'combined';
    task_runs   = 1:1:12;
    switch ptNum
        case 8010
            task_runs(task_runs==10) = [];
        case 8019
            task_runs(task_runs==10) = [];
    end
    
    fl_mot = dir(fullfile(func_dir, sprintf('s-%dmm_*task-RDMmotion*.nii', fwhm)))
    fl_col = dir(fullfile(func_dir, sprintf('s-%dmm_*task-RDMcolor*.nii', fwhm)))
    
    [mm,cc] = deal(1);
    clear fl
    for ii = 1:length(task_runs)
        if mod(task_runs(ii),2) == 1
            fl(ii) = fl_mot(mm);
            mm = mm+1;
        else
            fl(ii) = fl_col(cc);
            cc = cc+1;
        end
    end
    
    
    cl_mot = dir(sprintf('%s/sub*task-RDMmotion*.tsv', cfd_dir))
    cl_col = dir(sprintf('%s/sub*task-RDMcolor*.tsv', cfd_dir))
    
    [mm,cc] = deal(1);
    clear cl
    for ii = 1:length(task_runs)
        if mod(task_runs(ii),2) == 1
            cl(ii) = cl_mot(mm);
            mm = mm+1;
        else
            cl(ii) = cl_col(cc);
            cc = cc+1;
        end
    end
    
    
    
else % if COLOR
    
    movtName = 'color';
    task_runs   = 2:2:12;
    
    
    switch ptNum
        case 8010
            task_runs(task_runs==10) = [];
        case 8019
            task_runs(task_runs==10) = [];
    end
    
    cl = dir(sprintf('%s/sub*task-RDMcolor*.tsv', cfd_dir))
    fl = dir(fullfile(func_dir, sprintf('s-%dmm_*task-RDMcolor*.nii', fwhm)))
    
end





if exist(save_dir, 'dir')
    rmdir(save_dir, 's'); % remove existing results folder
end

% make directories
mkdir(save_dir); % save results
mkdir(cond_dir);  % save conditions




%% === load behavioral results
tl = dir(sprintf('%s/20*%d.mat', behav_dir, ptNum))
load(fullfile(tl(end).folder, tl(end).name));

r = data.results.rand.session;

performance = nanmean(r.behav.acc(:, task_runs))
task_runs(performance < .60) = []; %remove low performance runs


prevConf = [ones(1, size(r.confs,2)).*3; r.confs(1:end-1,:)];
nextConf = [r.confs(2:end,:); ones(1, size(r.confs,2)).*3];
taskMot    = 1-(sum(sum(isfinite(r.behav.rt(:, task_runs(mod(task_runs,2)==1))))) ./ sum(sum(isfinite(r.behav.rt))));
taskCol    = 1-(sum(sum(isfinite(r.behav.rt(:, task_runs(mod(task_runs,2)==0))))) ./ sum(sum(isfinite(r.behav.rt))));



if ~isempty(regexp(name, '8-2'))
    
    sd = load(fullfile(sense_dir, 'fit_senseDyn_8-2'))
    
elseif ~isempty(regexp(name, '16-4'))
    
    
    sd = load(fullfile(sense_dir, 'fit_senseDyn_16-4'))
    
elseif ~isempty(regexp(name, '32-8'))
    
    
    sd = load(fullfile(sense_dir, 'fit_senseDyn_32-8'))
    
    
elseif ~isempty(regexp(name, '10-60'))
    
    
    sd = load(fullfile(sense_dir, 'fit_senseFourier_10-60'))
    
    
elseif ~isempty(regexp(name, '10-100'))
    
    sd = load(fullfile(sense_dir, 'fit_senseFourier_10-100'))
    
end




%% === make condition files
if any(regName{2}) % if BLK
    
    
else % if MOTION or COLOR
end

nRuns = length(fl);

condList = [];
all_R    = [];
rlen     = [];

for rr = 1:nRuns
    
    
    
    %% === make task regressors
    center = @(x) x-nanmean(x);
    
    isRT    = isfinite(r.behav.rt(:, task_runs(rr)));
    acc     = r.behav.acc(:, task_runs(rr));
    
    %     distCenter = nanmean(r.confs(:,task_runs(rr)));
    %     targCenter = nanmean(r.attendCohs(:,task_runs(rr)));
    
    [durations, names, onsets, pmod, orth] = deal([]);
    
    
    
    trialDur = data.param.dotSettings.InterogationFixedTimeLength;
    switch confound
        
        case 'dur0'
            
            trialDur = 0;
            
    end
    
    
    
    switch name
        
            
        case 'cLinear'  % =====================================================
            
            
            runList = [1:6]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = -.5*(r.confs(isRT, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = .5*(r.attendCohs(isRT, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = (r.behav.resp(isRT, task_runs(rr))-2)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 5: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==5) = [];
            end
            
            nc = nc+1;
            
            
            % =============== 6: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==6) = [];
            end
            
            
        case 'cTargDist'  % =====================================================
            
            
            runList = [1:4]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(:, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = -.5*(r.confs(:, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = .5*(r.attendCohs(:, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: response
            resp = (r.behav.resp(:, task_runs(rr))-2)';
            resp(~isfinite(resp)) = 0;
            
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = resp;
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
        case 'cTargDistPerform'  % =====================================================
            
            
            runList = [1:7]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(:, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = -.5*(r.confs(:, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = .5*(r.attendCohs(:, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: response
            resp = (r.behav.resp(:, task_runs(rr))-2)';
            resp(~isfinite(resp)) = 0;
            
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = resp;
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 5: rt
            rt = center(r.behav.rt(:, task_runs(rr)))';
            rt(~isfinite(rt)) = 0;
            
            pmod(nc).name{pm}     = 'rt';
            pmod(nc).param{pm}    = rt;
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 6: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                acc = center(r.behav.acc(:, task_runs(rr)))';
                
                pmod(nc).name{pm}     = 'acc';
                pmod(nc).param{pm}    = acc;
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==6) = [];
            end
            
            
            
            % =============== 7: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                
                lapse = center(isnan(r.behav.rt(:, task_runs(rr)))');
                
                pmod(nc).name{pm}     = 'lapse';
                pmod(nc).param{pm}    = lapse;
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==7) = [];
            end
            
            
        case 'cTargDist-mot'  % =====================================================
            
            
            runList = [1:4]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(:, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = -.5*(r.confs(:, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = .5*(r.attendCohs(:, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: response
            resp = (r.behav.resp(:, task_runs(rr))-2)';
            resp(~isfinite(resp)) = 0;
            
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = resp;
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;

            
        case 'cLinear-mot'  % =====================================================
            
            
            runList = [1:6]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = -.5*(r.confs(isRT, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = .5*(r.attendCohs(isRT, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = (r.behav.resp(isRT, task_runs(rr))-2)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 5: error
            if ~isempty(r.time.dot_relStart(acc==0 & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==5) = [];
            end
            
            nc = nc+1;
            
            
            % =============== 6: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT & acc==0, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT & acc==0, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==6) = [];
            end
            
            
        case 'cLinear-comb'  % =====================================================
            
            
            runList = [1:6]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = -.5*(r.confs(isRT, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = .5*(r.attendCohs(isRT, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = (r.behav.resp(isRT, task_runs(rr))-2)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 5: error
            if ~isempty(r.time.dot_relStart(acc==0 & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==5) = [];
            end
            
            nc = nc+1;
            
            
            % =============== 6: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT & acc==0, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT & acc==0, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==6) = [];
            end
            
            
        case 'cLinearRT'  % =====================================================
            
            
            runList = [1:7]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = -.5*(r.confs(isRT, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = .5*(r.attendCohs(isRT, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = (r.behav.resp(isRT, task_runs(rr))-2)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 5: RT
            RT = center(r.behav.rt(isRT, task_runs(rr)))';
            
            pmod(nc).name{pm}     = 'RT';
            pmod(nc).param{pm}    = RT;
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 6: acc
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'acc';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==6) = [];
            end
            
            nc = nc+1;
            
            
            % =============== 7: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==7) = [];
            end
            
            
        case 'cLinearRTx'  % =====================================================
            
            
            runList = [1:9]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center(r.behav.resp(isRT, task_runs(rr))-2)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 5: RT
            RT = center(r.behav.rt(isRT, task_runs(rr)))';
            
            pmod(nc).name{pm}     = 'RT';
            pmod(nc).param{pm}    = RT;
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 6: dist x RT
            pmod(nc).name{pm}     = 'distRT';
            pmod(nc).param{pm}    = RT .* center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 7: targ x RT
            pmod(nc).name{pm}     = 'targRT';
            pmod(nc).param{pm}    = RT .* center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 8: acc
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'acc';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==8) = [];
            end
            
            nc = nc+1;
            
            
            % =============== 9: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==9) = [];
            end
            
            
        case 'cQuad'  % =====================================================
            
            dot_relStart = r.time.dot_start - r.time.blOnset;
            runList = [1:7]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials

            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: abs distractor congruence
            pmod(nc).name{pm}     = 'dist2';
            pmod(nc).param{pm}    = center((-.5*(r.confs(isRT, task_runs(rr)) - 3)').^2);
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 5: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center((r.behav.resp(isRT, task_runs(rr))-2)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 6: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==6) = [];
            end
            
            
            nc = nc+1;
            
            
            % =============== 7: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==7) = [];
            end
            
            
        case 'cQuadResp'  % =====================================================
            
            
            runList = [1:10]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'distR';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)') .* ((r.behav.resp(isRT, task_runs(rr))-2)' == 1);
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
             % 3: distractor congruence
            pmod(nc).name{pm}     = 'distL';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)') .* ((r.behav.resp(isRT, task_runs(rr))-2)' == -1);
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: abs distractor congruence
            pmod(nc).name{pm}     = 'dist2R';
            pmod(nc).param{pm}    = center((-.5*(r.confs(isRT, task_runs(rr)) - 3)').^2) .* ((r.behav.resp(isRT, task_runs(rr))-2)' == 1);
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
             % 5: abs distractor congruence
            pmod(nc).name{pm}     = 'dist2L';
            pmod(nc).param{pm}    = center((-.5*(r.confs(isRT, task_runs(rr)) - 3)').^2) .* ((r.behav.resp(isRT, task_runs(rr))-2)' == -1);
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 6: target coherence
            pmod(nc).name{pm}     = 'targetR';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)') .* ((r.behav.resp(isRT, task_runs(rr))-2)' == 1);
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 7: target coherence
            pmod(nc).name{pm}     = 'targetL';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)') .* ((r.behav.resp(isRT, task_runs(rr))-2)' == -1);
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 8: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center((r.behav.resp(isRT, task_runs(rr))-2)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            
            
            % 9: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==9) = [];
            end
            
            
            
            nc = nc+1;
            
            
            % =============== 10: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==10) = [];
            end

            
        case 'cSgnSal'  % =====================================================
            
            
            runList = [1:8]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: abs distractor congruence
            pmod(nc).name{pm}     = 'dsal';
            pmod(nc).param{pm}    = center(abs(-.5*(r.confs(isRT, task_runs(rr)) - 3)'));
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: sgn distractor congruence
            pmod(nc).name{pm}     = 'dsgn';
            pmod(nc).param{pm}    = center(sign(-.5*(r.confs(isRT, task_runs(rr)) - 3)'));
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 5: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 6: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center((r.behav.resp(isRT, task_runs(rr))-2)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 7: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==7) = [];
            end
            
            
            
            nc = nc+1;
            
            
            % =============== 8: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==8) = [];
            end
            
            
        case 'cHiLo'  % =====================================================
            
            
            runList = [1:7]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            
            % 2: distractor congruence HI
            pmod(nc).name{pm}     = 'distHi';
            pmod(nc).param{pm}    = -sign(ismember(r.confs(isRT, task_runs(rr)), [1,5]) .* (r.confs(isRT, task_runs(rr)) - 3))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: distractor congruence LO
            pmod(nc).name{pm}     = 'distLo';
            pmod(nc).param{pm}    = -sign(ismember(r.confs(isRT, task_runs(rr)), [2,4]) .* (r.confs(isRT, task_runs(rr)) - 3))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = .5*(r.attendCohs(isRT, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 5: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = (r.behav.resp(isRT, task_runs(rr))-2)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 6: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==6) = [];
            end
            
            nc = nc+1;
            
            
            % =============== 7: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==7) = [];
            end
            
            
        case 'cHiLoSal'  % =====================================================
            
            
            runList = [1:9]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            
            % 2: distractor congruence HI
            pmod(nc).name{pm}     = 'distHi';
            pmod(nc).param{pm}    = -sign(ismember(r.confs(isRT, task_runs(rr)), [1,5]) .* (r.confs(isRT, task_runs(rr)) - 3))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: distractor congruence LO
            pmod(nc).name{pm}     = 'distLo';
            pmod(nc).param{pm}    = -sign(ismember(r.confs(isRT, task_runs(rr)), [2,4]) .* (r.confs(isRT, task_runs(rr)) - 3))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: distractor Salience HI
            pmod(nc).name{pm}     = 'distSalHi';
            pmod(nc).param{pm}    = ismember(r.confs(isRT, task_runs(rr)), [1,5])';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 5: distractor Salience LO
            pmod(nc).name{pm}     = 'distSalLo';
            pmod(nc).param{pm}    =  ismember(r.confs(isRT, task_runs(rr)), [2,4])';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 6: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = .5*(r.attendCohs(isRT, task_runs(rr)) - 3)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 7: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = (r.behav.resp(isRT, task_runs(rr))-2)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 8: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==8) = [];
            end
            
            nc = nc+1;
            
            
            % =============== 9: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==9) = [];
            end
            
            
        case 'cAdapt'   % =====================================================
            
            
            
            runList = [1:8]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            pm = 1;
            
            
            % 2: distractor congruence (accurate trials)
            pmod(nc).name{pm}     = 'dist0';
            pmod(nc).param{pm}    = center(-.5*((r.confs(isRT, task_runs(rr)) - 3))');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: distractor congruence (accurate trials)
            pmod(nc).name{pm}     = 'dist1';
            pmod(nc).param{pm}    = center(-.5*((prevConf(isRT, task_runs(rr)) - 3))');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: distractor congruence (accurate trials)
            pmod(nc).name{pm}     = 'adapt';
            pmod(nc).param{pm}    = center(-.5*((r.confs(isRT, task_runs(rr)) - 3))') .* center(-.5*((prevConf(isRT, task_runs(rr)) - 3))');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 5: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 6: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center(r.behav.resp(isRT, task_runs(rr))-2)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 7: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==7) = [];
            end
            
            nc = nc+1;
            
            
            % =============== 8: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==8) = [];
            end
            
            
        case 'cAdaptContrast'   % =====================================================
            
            
            
            runList = [1:10]; % keep track of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            pm = 1;
            
            
            % 2: distractor congruence (accurate trials)
            pmod(nc).name{pm}     = 'dist0';
            pmod(nc).param{pm}    = center(-.5*((r.confs(isRT, task_runs(rr)) - 3))');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: distractor congruence (accurate trials)
            pmod(nc).name{pm}     = 'distPast';
            pmod(nc).param{pm}    = center(-.5*((prevConf(isRT, task_runs(rr)) - 3))');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: distractor congruence (accurate trials)
            pmod(nc).name{pm}     = 'distFuture';
            pmod(nc).param{pm}    = center(-.5*((nextConf(isRT, task_runs(rr)) - 3))');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            
            % 5: distractor congruence (accurate trials)
            pmod(nc).name{pm}     = 'adaptPast';
            pmod(nc).param{pm}    = center(-.5*((r.confs(isRT, task_runs(rr)) - 3))') .* center(-.5*((prevConf(isRT, task_runs(rr)) - 3))');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
             % 6: distractor congruence (accurate trials)
            pmod(nc).name{pm}     = 'adaptFuture';
            pmod(nc).param{pm}    = center(-.5*((r.confs(isRT, task_runs(rr)) - 3))') .* center(-.5*((nextConf(isRT, task_runs(rr)) - 3))');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 7: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 8: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center(r.behav.resp(isRT, task_runs(rr))-2)';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 9: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==9) = [];
            end
            
            
            nc = nc+1;
            
            
            % =============== 10: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==10) = [];
            end

            
        case 'cSenseFourier_10-60'  % =====================================================
            
            ptSel = find(sd.pts == ptNum);
            distGain    = sd.out(ptSel).distGain_rt(:, sd.out(ptSel).colRun == task_runs(rr));
            targGain    = sd.out(ptSel).targGain_rt(:, sd.out(ptSel).colRun == task_runs(rr));
%             distVel     = sd.out(ptSel).distVel_rt(:, sd.out(ptSel).colRun == task_runs(rr));
%             targVel     = sd.out(ptSel).targVel_rt(:, sd.out(ptSel).colRun == task_runs(rr));
            
            
            runList = [1:16]; % keep track of number of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            pm = 1;
            
            
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            
            
            
            % 3: distractor gain 10
            pmod(nc).name{pm}     = 'distGain10';
            pmod(nc).param{pm}    = center(distGain(1, isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: distractor gain 20
            pmod(nc).name{pm}     = 'distGain20';
            pmod(nc).param{pm}    = center(distGain(2, isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 5: distractor gain 30
            pmod(nc).name{pm}     = 'distGain30';
            pmod(nc).param{pm}    = center(distGain(3, isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 6: distractor gain 40
            pmod(nc).name{pm}     = 'distGain40';
            pmod(nc).param{pm}    = center(distGain(4, isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 7: distractor gain 50
            pmod(nc).name{pm}     = 'distGain50';
            pmod(nc).param{pm}    = center(distGain(5, isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            
            
            
            
            
            
            
            % 8: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 9: distractor gain 10
            pmod(nc).name{pm}     = 'targGain10';
            pmod(nc).param{pm}    = center(targGain(1, isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 10: distractor gain 10
            pmod(nc).name{pm}     = 'targGain20';
            pmod(nc).param{pm}    = center(targGain(2, isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 11: distractor gain 10
            pmod(nc).name{pm}     = 'targGain30';
            pmod(nc).param{pm}    = center(targGain(3, isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 12: distractor gain 10
            pmod(nc).name{pm}     = 'targGain40';
            pmod(nc).param{pm}    = center(targGain(4, isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 13: distractor gain 10
            pmod(nc).name{pm}     = 'targGain50';
            pmod(nc).param{pm}    = center(targGain(5, isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            
            % 14: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center((r.behav.resp(isRT, task_runs(rr))-2)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 15: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==15) = [];
            end
            
            
            
            nc = nc+1;
            
            
            % =============== 16: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==16) = [];
            end

            
        case 'cSenseFourierAll_10-100'  % =====================================================
            
            
            ptSel = find(sd.pts == ptNum);
            
            distGain_rt    = sd.out(ptSel).distGainAll_rt(sd.out(ptSel).colRun == task_runs(rr));
            distVel_rt     = sd.out(ptSel).distVelAll_rt(sd.out(ptSel).colRun == task_runs(rr));
            distGain_acc    = sd.out(ptSel).distGainAll_acc(sd.out(ptSel).colRun == task_runs(rr));
            distVel_acc     = sd.out(ptSel).distVelAll_acc(sd.out(ptSel).colRun == task_runs(rr));
            
            targGain_rt    = sd.out(ptSel).targGainAll_rt(sd.out(ptSel).colRun == task_runs(rr));
            targVel_rt     = sd.out(ptSel).targVelAll_rt(sd.out(ptSel).colRun == task_runs(rr));
            targGain_acc    = sd.out(ptSel).targGainAll_acc(sd.out(ptSel).colRun == task_runs(rr));
            targVel_acc     = sd.out(ptSel).targVelAll_acc(sd.out(ptSel).colRun == task_runs(rr));
            
            runList = [1:15]; % keep track of number of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            pm = 1;
            
            
            
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: distractor gain
            pmod(nc).name{pm}     = 'distGain_{rt}';
            pmod(nc).param{pm}    = center(distGain_rt(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: distractor gain Vel
            pmod(nc).name{pm}     = 'distVel_{rt}';
            pmod(nc).param{pm}    = center(distVel_rt(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 5: distractor gain
            pmod(nc).name{pm}     = 'distGain_{acc}';
            pmod(nc).param{pm}    = center(distGain_acc(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 6: distractor gain Vel
            pmod(nc).name{pm}     = 'distVel_{acc}';
            pmod(nc).param{pm}    = center(distVel_acc(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            
            
            
            
            % 7: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 8: distractor gain
            pmod(nc).name{pm}     = 'targGain_{rt}';
            pmod(nc).param{pm}    = center(targGain_rt(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 9: distractor gain vel
            pmod(nc).name{pm}     = 'targVel_{rt}';
            pmod(nc).param{pm}    = center(targVel_rt(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 10: distractor gain
            pmod(nc).name{pm}     = 'targGain_{acc}';
            pmod(nc).param{pm}    = center(targGain_acc(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 11: distractor gain vel
            pmod(nc).name{pm}     = 'targVel_{acc}';
            pmod(nc).param{pm}    = center(targVel_acc(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            
            
            
            % 12: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center((r.behav.resp(isRT, task_runs(rr))-2)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 13: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==13) = [];
            end
            
            
            % 14: RT
            pmod(nc).name{pm}     = 'RT';
            pmod(nc).param{pm}    = center(r.behav.rt(isRT, task_runs(rr)))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            
            % =============== 15: lapse trials
            nc = nc+1;
            
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==15) = [];
            end
            
              
        case 'cSenseDyn_8-2'  % =====================================================
            
            ptSel = find(sd.pts == ptNum);
            distGain    = sd.out(ptSel).distGain_rt(sd.out(ptSel).colRun == task_runs(rr));
            distVel     = sd.out(ptSel).distVel_rt(sd.out(ptSel).colRun == task_runs(rr));
            targGain    = sd.out(ptSel).targGain_rt(sd.out(ptSel).colRun == task_runs(rr));
            targVel     = sd.out(ptSel).targVel_rt(sd.out(ptSel).colRun == task_runs(rr));
            
            
            runList = [1:10]; % keep track of number of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: distractor gain
            pmod(nc).name{pm}     = 'distGain';
            pmod(nc).param{pm}    = center(distGain(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: distractor vel
            pmod(nc).name{pm}     = 'distVel';
            pmod(nc).param{pm}    = center(distVel(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 5: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 6: target gain
            pmod(nc).name{pm}     = 'targGain';
            pmod(nc).param{pm}    = center(targGain(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 7: target vel
            pmod(nc).name{pm}     = 'targVel';
            pmod(nc).param{pm}    = center(targVel(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 8: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center((r.behav.resp(isRT, task_runs(rr))-2)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 9: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==9) = [];
            end
            
            
            
            nc = nc+1;
            
            
            % =============== 10: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==10) = [];
            end
            
            
        case 'cSenseDyn_16-4'  % =====================================================
            
            ptSel = find(sd.pts == ptNum);
            distGain    = sd.out(ptSel).distGain_rt(sd.out(ptSel).colRun == task_runs(rr));
            distVel     = sd.out(ptSel).distVel_rt(sd.out(ptSel).colRun == task_runs(rr));
            targGain    = sd.out(ptSel).targGain_rt(sd.out(ptSel).colRun == task_runs(rr));
            targVel     = sd.out(ptSel).targVel_rt(sd.out(ptSel).colRun == task_runs(rr));
            
            
            runList = [1:10]; % keep track of number of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: distractor gain
            pmod(nc).name{pm}     = 'distGain';
            pmod(nc).param{pm}    = center(distGain(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: distractor vel
            pmod(nc).name{pm}     = 'distVel';
            pmod(nc).param{pm}    = center(distVel(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 5: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 6: target gain
            pmod(nc).name{pm}     = 'targGain';
            pmod(nc).param{pm}    = center(targGain(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 7: target vel
            pmod(nc).name{pm}     = 'targVel';
            pmod(nc).param{pm}    = center(targVel(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 8: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center((r.behav.resp(isRT, task_runs(rr))-2)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 9: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==9) = [];
            end
            
            
            
            nc = nc+1;
            
            
            % =============== 10: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==10) = [];
            end
            
            
        case 'cSenseDyn_32-8'  % =====================================================
            
            ptSel = find(sd.pts == ptNum);
            distGain    = sd.out(ptSel).distGain_rt(sd.out(ptSel).colRun == task_runs(rr));
            distVel     = sd.out(ptSel).distVel_rt(sd.out(ptSel).colRun == task_runs(rr));
            targGain    = sd.out(ptSel).targGain_rt(sd.out(ptSel).colRun == task_runs(rr));
            targVel     = sd.out(ptSel).targVel_rt(sd.out(ptSel).colRun == task_runs(rr));
            
            
            runList = [1:10]; % keep track of number of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: distractor gain
            pmod(nc).name{pm}     = 'distGain';
            pmod(nc).param{pm}    = center(distGain(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: distractor vel
            pmod(nc).name{pm}     = 'distVel';
            pmod(nc).param{pm}    = center(distVel(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 5: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 6: target gain
            pmod(nc).name{pm}     = 'targGain';
            pmod(nc).param{pm}    = center(targGain(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 7: target vel
            pmod(nc).name{pm}     = 'targVel';
            pmod(nc).param{pm}    = center(targVel(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 8: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center((r.behav.resp(isRT, task_runs(rr))-2)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 9: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==9) = [];
            end
            
            
            
            nc = nc+1;
            
            
            % =============== 10: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==10) = [];
            end
            
            
        case 'cSenseDyn_X_16-4'  % =====================================================
            
            ptSel = find(sd.pts == ptNum);
            distGain    = sd.out(ptSel).distGain_rt(sd.out(ptSel).colRun == task_runs(rr));
            distVel     = sd.out(ptSel).distVel_rt(sd.out(ptSel).colRun == task_runs(rr));
            targGain    = sd.out(ptSel).targGain_rt(sd.out(ptSel).colRun == task_runs(rr));
            targVel     = sd.out(ptSel).targVel_rt(sd.out(ptSel).colRun == task_runs(rr));
            
            
            runList = [1:11]; % keep track of number of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: distractor gain
            pmod(nc).name{pm}     = 'distGain';
            pmod(nc).param{pm}    = center(distGain(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: distractor X
            pmod(nc).name{pm}     = 'distX';
            pmod(nc).param{pm}    = center(distGain(isRT))' .* center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 5: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 6: target gain
            pmod(nc).name{pm}     = 'targGain';
            pmod(nc).param{pm}    = center(targGain(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 7: target X
            pmod(nc).name{pm}     = 'targX';
            pmod(nc).param{pm}    = center(targGain(isRT))' .* center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 8: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center((r.behav.resp(isRT, task_runs(rr))-2)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 9: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==9) = [];
            end
            
            
            
            nc = nc+1;
            
            
            % =============== 10: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==10) = [];
            end
            
            
        case 'cSenseDyn_RTACC_16-4'  % =====================================================
            
            ptSel = find(sd.pts == ptNum);
            distGain_rt    = sd.out(ptSel).distGain_rt(sd.out(ptSel).colRun == task_runs(rr));
            targGain_rt    = sd.out(ptSel).targGain_rt(sd.out(ptSel).colRun == task_runs(rr));
            distGain_acc    = sd.out(ptSel).distGain_acc(sd.out(ptSel).colRun == task_runs(rr));
            targGain_acc    = sd.out(ptSel).targGain_acc(sd.out(ptSel).colRun == task_runs(rr));
            
            
            runList = [1:11]; % keep track of number of conditions
            
            
            nc = 1;
            
            %  =============== 1: response trials
            durations{nc}        = trialDur;
            names{nc}            = 'trial';
            onsets{nc}           = r.time.dot_relStart(isRT, task_runs(rr))';
            orth{nc}             = 0;
            
            pm = 1;
            
            % 2: distractor congruence
            pmod(nc).name{pm}     = 'dist';
            pmod(nc).param{pm}    = center(-.5*(r.confs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 3: distractor gain
            pmod(nc).name{pm}     = 'distGainRT';
            pmod(nc).param{pm}    = center(distGain_rt(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 4: distractor X
            pmod(nc).name{pm}     = 'distGainAcc';
            pmod(nc).param{pm}    = -center(distGain_acc(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 5: target coherence
            pmod(nc).name{pm}     = 'target';
            pmod(nc).param{pm}    = center(.5*(r.attendCohs(isRT, task_runs(rr)) - 3)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 6: target gain
            pmod(nc).name{pm}     = 'targGainRT';
            pmod(nc).param{pm}    = center(targGain_rt(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 7: target X
            pmod(nc).name{pm}     = 'targGainAcc';
            pmod(nc).param{pm}    = -center(targGain_acc(isRT))';
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            % 8: response
            pmod(nc).name{pm}     = 'resp';
            pmod(nc).param{pm}    = center((r.behav.resp(isRT, task_runs(rr))-2)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;
            
            
            
            % 9: error
            if ~isempty(r.time.dot_relStart(~acc & isRT, task_runs(rr))')
                
                pmod(nc).name{pm}     = 'error';
                pmod(nc).param{pm}    = center(r.behav.acc(isRT, task_runs(rr)))';
                pmod(nc).poly{pm}     = 1;
                pm = pm+1;
                
            else
                runList(runList==9) = [];
            end
            
            
            
            nc = nc+1;
            
            
            % =============== 10: lapse trials
            if ~isempty(r.time.dot_relStart(~isRT, task_runs(rr))')
                durations{nc}        = trialDur;
                names{nc}            = 'lapse';
                onsets{nc}           = r.time.dot_relStart(~isRT, task_runs(rr))';
                orth{nc}             = 0;
                
                nc = nc+1;
            else
                runList(runList==10) = [];
            end

            
    end
    
    
    
    
    % save mat
    save(sprintf('%s/%s_run%d', cond_dir, analysis, rr),...
        'names', 'onsets', 'durations', 'pmod', 'orth')
    
    
    
    %% === make confound regressors
    
    
    cs = tdfread(fullfile(cl(rr).folder, cl(rr).name));
    
    % framewise displacement
    fd = [0; str2num(cs.framewise_displacement(2:end,:))];
    
    % movement deriv
    trans_x_derivative1 = [0; str2num(cs.trans_x_derivative1(2:end,:))];
    trans_y_derivative1 = [0; str2num(cs.trans_y_derivative1(2:end,:))];
    trans_z_derivative1 = [0; str2num(cs.trans_z_derivative1(2:end,:))];
    rot_x_derivative1 = [0; str2num(cs.rot_x_derivative1(2:end,:))];
    rot_y_derivative1 = [0; str2num(cs.rot_y_derivative1(2:end,:))];
    rot_z_derivative1 = [0; str2num(cs.rot_z_derivative1(2:end,:))];
    
    trans_x_derivative1_power2 = [0; str2num(cs.trans_x_derivative1_power2(2:end,:))];
    trans_y_derivative1_power2 = [0; str2num(cs.trans_y_derivative1_power2(2:end,:))];
    trans_z_derivative1_power2 = [0; str2num(cs.trans_z_derivative1_power2(2:end,:))];
    rot_x_derivative1_power2 = [0; str2num(cs.rot_x_derivative1_power2(2:end,:))];
    rot_y_derivative1_power2 = [0; str2num(cs.rot_y_derivative1_power2(2:end,:))];
    rot_z_derivative1_power2 = [0; str2num(cs.rot_z_derivative1_power2(2:end,:))];
    
    
    R = ([...
        cs.trans_x, cs.trans_y, cs.trans_z,...
        cs.rot_x, cs.rot_y, cs.rot_z,...
        ...
        trans_x_derivative1, trans_y_derivative1, trans_z_derivative1,...
        rot_x_derivative1, rot_y_derivative1, rot_z_derivative1,...
        ...
        cs.trans_x_power2, cs.trans_y_power2, cs.trans_z_power2,...
        cs.rot_x_power2, cs.rot_y_power2, cs.rot_z_power2,...
        ...
        trans_x_derivative1_power2, trans_y_derivative1_power2, trans_z_derivative1_power2,...
        rot_x_derivative1_power2, rot_y_derivative1_power2, rot_z_derivative1_power2,...
        ]);
    
    
    names = { ...
        'tx', 'ty', 'tz',...
        'rx', 'ry', 'rz',...
        ...
        'dtx', 'dty', 'dtz',...
        'drx', 'dry', 'drz',...
        ...
        'tx2', 'ty2', 'tz2',...
        'rx2', 'ry2', 'rz2',...
        ...
        'dtx2', 'dty2', 'dtz2',...
        'drx2', 'dry2', 'drz2',...
        };



% cs = tdfread(fullfile(cl(rr).folder, cl(rr).name));
% 
% % framewise displacement
% fd = [0; str2num(cs.framewise_displacement(2:end,:))];
% 
% % movement deriv
% trans_x_derivative1 = [0; str2num(cs.trans_x_derivative1(2:end,:))];
% trans_y_derivative1 = [0; str2num(cs.trans_y_derivative1(2:end,:))];
% trans_z_derivative1 = [0; str2num(cs.trans_z_derivative1(2:end,:))];
% rot_x_derivative1 = [0; str2num(cs.rot_x_derivative1(2:end,:))];
% rot_y_derivative1 = [0; str2num(cs.rot_y_derivative1(2:end,:))];
% rot_z_derivative1 = [0; str2num(cs.rot_z_derivative1(2:end,:))];
% 
% trans_x_derivative1_power2 = [0; str2num(cs.trans_x_derivative1_power2(2:end,:))];
% trans_y_derivative1_power2 = [0; str2num(cs.trans_y_derivative1_power2(2:end,:))];
% trans_z_derivative1_power2 = [0; str2num(cs.trans_z_derivative1_power2(2:end,:))];
% rot_x_derivative1_power2 = [0; str2num(cs.rot_x_derivative1_power2(2:end,:))];
% rot_y_derivative1_power2 = [0; str2num(cs.rot_y_derivative1_power2(2:end,:))];
% rot_z_derivative1_power2 = [0; str2num(cs.rot_z_derivative1_power2(2:end,:))];
% 
% 
% R = ([...
%     fd,...
%     cs.trans_x, cs.trans_y, cs.trans_z,...
%     cs.rot_x, cs.rot_y, cs.rot_z,...
%     ]);
% 
% names = { ...
%     'fd',...
%     'tx', 'ty', 'tz',...
%     'rx', 'ry', 'rz',...
%     };
% 
% 
% find_fd = find(fd>1)
% censor = [];
% 
% if ~isempty(find_fd)
%     
%     for ii = 1:length(find_fd)
%         
%         censor = blkdiag(censor, fd == fd(find_fd(ii)));
%         names{end+1} = sprintf('cens%d',ii);
%         
%     end
%     
%     R = [R, censor];
%     
%     
% end


















% save & update list
save(sprintf('%s/confound_run%d', cfd_dir, rr),...
    'R', 'names');


runList = [runList, zeros(size(names))]
condList = [condList, runList];

all_R = [all_R; R];

rlen(rr) = size(R,1);




    
    
end

save(sprintf('%s/%s_condList', cond_dir, analysis), 'condList', 'task_runs', 'taskMot', 'taskCol');



%% plot motion
mkdir(sprintf('%s/spm-data/motion', root_dir));

f=figure; hold on;
f.Position = [100 100 1300 1000];

Rs = cumsum(rlen(1:end-1));


% % % translation % % % 
subplot(3,1,1); hold on;
plot(all_R(:,ismember(names, {'tx', 'ty', 'tz'})), '-', 'LineWidth', 2);
for rr = 1:length(Rs)
    xline(Rs(rr), '-k', 'LineWidth', .75);
end
legend({'tx', 'ty', 'tz'}, 'Location', 'northeastoutside')
title(sprintf('translation - sub-%d', ptNum))
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlabel('TR')


% % % rotation % % % 
subplot(3,1,2); hold on;
plot(all_R(:,ismember(names, {'rx', 'ry', 'rz'})), '-', 'LineWidth', 2);
for rr = 1:length(Rs)
    xline(Rs(rr), '-k', 'LineWidth', .75);
end
legend({'rx', 'ry', 'rz'}, 'Location', 'northeastoutside')
title(sprintf('rotation - sub-%d', ptNum))
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlabel('TR')


% % % % displacement % % % 
% subplot(3,1,3); hold on;
% plot(all_R(:,ismember(names, {'fd'})), '-', 'LineWidth', 2);
% for rr = 1:length(Rs)
%     xline(Rs(rr), '-k', 'LineWidth', .75);
% end
% legend({'fd'}, 'Location', 'northeastoutside')
% title(sprintf('displacement - sub-%d', ptNum))
% set(gca, 'TickDir', 'out', 'LineWidth', 1);
% xlabel('TR')


% % % wm/csf % % % 
% subplot(4,1,4); hold on;
% plot(all_R(:,ismember(names, {'csfwm'})), '-', 'LineWidth', 2);
% for rr = 1:length(Rs)
%     xline(Rs(rr), '-k', 'LineWidth', .75);
% end
% legend({'csfwm'}, 'Location', 'northeastoutside')
% title(sprintf('translation - sub-%d', ptNum))
% set(gca, 'TickDir', 'out', 'LineWidth', 1);
% xlabel('TR')




% save & close
saveas(f,sprintf('%s/spm-data/motion/%s_sub-%d.png', root_dir, movtName, ptNum))

close(f);


%% === get mask

ml = dir(sprintf('%s/anat/%s', pt_dir, mask_filt));
mask_file = fullfile(ml(end).folder, ml(end).name)



%% === BATCH ====
%  ==============

spm_jobman('initcfg')

% specify model for each run
for rr = 1:nRuns
    
    % get frames
    runFile = fullfile(fl(rr).folder, fl(rr).name)
    matlabbatch{rr}.spm.util.exp_frames.files = {runFile};
    matlabbatch{rr}.spm.util.exp_frames.frames = Inf;
    
    
    % add scan/cond info to run-specific session
    matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess(rr).scans(1) = ...
        cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{rr}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    
    matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess(rr).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    
    matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess(rr).multi = ...
        {sprintf('%s/%s_run%d.mat', cond_dir,analysis, rr)};
    
    matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess(rr).regress = struct('name', {}, 'val', {});
    matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess(rr).multi_reg = {sprintf('%s/confound_run%d.mat', cfd_dir, rr)};
    matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess(rr).hpf = 128;
    
end



% specify general model features
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.dir = {save_dir};
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.timing.units = 'secs';
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.timing.RT = data.param.TR;
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.timing.fmri_t = 16;
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.timing.fmri_t0 = 8;

matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{nRuns+1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.volt = 1;
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.global = 'None';
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.mthresh = 0.8;
% matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.mask = {mask_file};
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.cvi = cvi;



% estimate model
matlabbatch{nRuns+2}.spm.tools.rwls.fmri_rwls_est.spmmat(1) = ...
    cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{nRuns+1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{nRuns+2}.spm.stats.fmri_est.write_residuals = 1;
matlabbatch{nRuns+2}.spm.tools.rwls.fmri_rwls_est.method.Classical = 1;



%% === RUN

tic
spm_jobman('run', matlabbatch);
toc


end