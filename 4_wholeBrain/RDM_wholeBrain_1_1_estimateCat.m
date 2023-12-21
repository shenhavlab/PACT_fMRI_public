function RDM_wholeBrain_1_1_estimateCat(root_dir, spm_dir, ptNum, varargin)
%% first level analysis ===================================================
% Harrison Ritz 2021



%% === initialize
addpath(genpath(spm_dir)); % add spm12 folder
spm('defaults', 'fmri')


% set default RAM & use RAM for analysis
spm_get_defaults('maxmem',32 * 2^30)
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

behav_dir   = sprintf('%s/behavior', root_dir)
pt_dir      = sprintf('%s/spm-data/sub-%d', root_dir, ptNum)
func_dir    = sprintf('%s/func', pt_dir)
save_dir    = sprintf('%s/level 1/%s', pt_dir, analysis)
cond_dir    = sprintf('%s/conditions', pt_dir)
cfd_dir     = sprintf('%s/confound', pt_dir)
sense_dir   = fullfile(root_dir, 'RDM_fmri_scripts', 'behavior')

mask_filt   = 'sub*desc-brain_mask.nii'   % mask filter
grey_filt   = 'sub*label-GM_probseg.nii'   % mask filter


regName = regexp(name, {'mot', 'comb', 'localizer'});

if any(regName{1}) % if Motion

    movtName = 'motion';
    task_runs   = 1:2:12;


    switch ptNum
        case 8010
            task_runs(task_runs==10) = [];
        case 8019
            task_runs(task_runs==10) = [];
    end

    blkDur = 125*1.2;
    blkVec = cumsum([0, repmat(blkDur, [1, length(task_runs)-1])]);

    cl = dir(sprintf('%s/sub*task-RDMmotion*.tsv', cfd_dir))
    fl = dir(fullfile(func_dir, sprintf('s-%dmm_*task-RDMmotion*.nii', fwhm)))

elseif any(regName{2}) % if COMBINED

    movtName = 'combined';
    task_runs   = 1:1:12;

    colSel = 2:2:12;
    motSel = 1:2:12;

    switch ptNum
        case 8010
            task_runs(task_runs==10) = [];
            colSel(colSel == 10) = [];
            task_runs(task_runs==9) = [];
            motSel(motSel == 9) = [];
        case 8019
            task_runs(task_runs==10) = [];
            colSel(colSel == 10) = [];
            task_runs(task_runs==9) = [];
            motSel(motSel == 9) = [];
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

    blkVec = [0, 1:length(task_runs)-1];
    blkVec(ismember(blkVec, colSel)) = 392*1.2;
    blkVec(ismember(blkVec, motSel)) = 125*1.2;
    blkVec = cumsum(blkVec);



elseif any(regName{3}) % if localizer

    movtName = 'localizer';
    task_runs   = 1:2;

    loc_ords = [1,1,0,0,1,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0];

    switch ptNum
        case 8006
            task_runs = 1;
    end

    blkDur = 195*1.2;
    blkVec = cumsum([0, repmat(blkDur, [1, length(task_runs)-1])]);

    cl = dir(sprintf('%s/sub*task-RDMlocalizer*.tsv', cfd_dir))
    fl = dir(fullfile(func_dir, sprintf('s-%dmm_*task-RDMlocalizer*.nii', fwhm)))



else % if COLOR

    movtName = 'color';
    task_runs   = 2:2:12;

    switch ptNum
        case 8010
            task_runs(task_runs==10) = [];
        case 8019
            task_runs(task_runs==10) = [];
    end

    blkDur = 392*1.2;
    blkVec = cumsum([0, repmat(blkDur, [1, length(task_runs)-1])]);


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

% filter by performance
% performance = nanmean(r.behav.acc(:, task_runs))
% if any(performance < .60)
%     disp(' ======= REMOVING BAD BLOCK ========= ')
%     disp(performance)
%     disp(' ======= ================== ========= ')
% end
%
% task_runs(performance < .60) = []; %remove low performance runs


taskMot    = 1-(sum(sum(isfinite(r.behav.rt(:, task_runs(mod(task_runs,2)==1))))) ./ sum(sum(isfinite(r.behav.rt))));
taskCol    = 1-(sum(sum(isfinite(r.behav.rt(:, task_runs(mod(task_runs,2)==0))))) ./ sum(sum(isfinite(r.behav.rt))));

prevDist = [ones(1, size(r.confs,2)).*3; r.confs(1:end-1,:)];
nextDist = [r.confs(2:end,:); ones(1, size(r.confs,2)).*3];

prev2Dist = [ones(2, size(r.confs,2)).*3; r.confs(1:end-2,:)];
next2Dist = [r.confs(3:end,:); ones(2, size(r.confs,2)).*3];

prevTarg = [ones(1, size(r.attendCohs,2)).*3; r.attendCohs(1:end-1,:)];
nextTarg = [r.attendCohs(2:end,:); ones(1, size(r.attendCohs,2)).*3];

prev2Targ = [ones(2, size(r.attendCohs,2)).*3; r.confs(1:end-2,:)];
next2Targ = [r.attendCohs(3:end,:); ones(2, size(r.confs,2)).*3];


prevRT = [ones(1, size(r.confs,2)).*nanmean(r.behav.rt(1:end-1,:)); r.behav.rt(1:end-1,:)];
nextRT = [r.behav.rt(2:end,:); ones(1, size(r.confs,2)).*nanmean(r.behav.rt(2:end,:))];

prevAcc = [ones(1, size(r.confs,2)).*nanmean(r.behav.acc(1:end-1,:));   r.behav.acc(1:end-1,:)];
nextAcc = [r.behav.acc(2:end,:);   ones(1, size(r.confs,2)).*nanmean(r.behav.acc(2:end,:))];



prevResp = [ones(1, size(r.confs,2)).*2;   r.behav.resp(1:end-1,:)];
nextResp = [r.behav.resp(2:end,:);   ones(1, size(r.confs,2)).*2];




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

% for rr = 1:nRuns



%% === make task regressors

vec = @(x) x(:);
center = @(x) x-nanmean(x);
nanz = @(x) (x-nanmean(x))./nanstd(x);
selz = @(x, sel) ((x-nanmean(x(logical(sel)))) ./ nanstd(x(logical(sel)))).*sel;
selc = @(x, sel) (x-nanmean(x(logical(sel)))).*sel;

[durations, names, onsets, pmod, orth, runList] = deal([]);

trialDur = data.param.dotSettings.InterogationFixedTimeLength;
dot_allStart = r.time.dot_relStart(:, task_runs) + blkVec;


switch confound

    case 'dur0'

        trialDur = 0;

end



switch name


    %% ===   JUNE 2022   ============================================================================================



    case 'featureBlk'    % =====================================================
        %%


        %  ---- resp -----
        resp = r.behav.resp(:, task_runs);
        resp(resp == 1) = -1;
        resp(resp == 3) = 1;

        corResp = r.task.corrResp(:, task_runs);
        corResp(corResp == 1) = -1;
        corResp(corResp == 3) = 1;

        %  ---- RT -----
        rt = r.behav.rt(:, task_runs);
        prevRT = prevRT(:, task_runs);
        nextRT = nextRT(:, task_runs);

        %  ---- Acc -----
        acc = r.behav.acc(:, task_runs);
        prevAcc = prevAcc(:, task_runs);
        nextAcc = nextAcc(:, task_runs);


        %  ---- Dist -----
        dist = r.confs(:, task_runs);
        prevDist = prevDist(:, task_runs);
        nextDist = nextDist(:, task_runs);
        prev2Dist = prev2Dist(:, task_runs);
        next2Dist = next2Dist(:, task_runs);

        %  ---- Targ -----
        targ = r.attendCohs(:, task_runs);
        prevTarg = prevTarg(:, task_runs);
        nextTarg = nextTarg(:, task_runs);
        prev2Targ = prev2Targ(:, task_runs);
        next2Targ = next2Targ(:, task_runs);




        sel = acc==1;


        nc = 1;

        %  =============== 1: response trials

        durations{nc}        = trialDur;
        names{nc}            = 'trial';
        onsets{nc}           = dot_allStart(sel)';
        orth{nc}             = 0;

        pm = 1;



        for bb = 1:size(acc,2)


            blksel = deal(zeros(size(acc)));
            blksel(:, bb) = 1;


            % targ resp
            pmod(nc).name{pm}     = sprintf('targResp_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist resp
            pmod(nc).name{pm}     = sprintf('distResp_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targCoh_%d', bb);
            pmod(nc).param{pm}    = selz((targ(sel)-3)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('distCoh_%d', bb);
            pmod(nc).param{pm}    = selz(abs(3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('cong_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;



        end




        runList = [1:(length(pmod(1).name)+2)] % keep track of conditions

        % ===== COLINEARITY =====
        xs = []; hx = [];
        for ii = 1:(pm-1)
            xs(:,ii) = pmod(nc).param{ii}';
            hx(:,ii) = conv(pmod(nc).param{ii}', spm_hrf(1.5));
        end

        disp(xs(1:50,1:5))
        statX   = [mean(hx); std(hx)]
        %         corrX   = corr(hx)
        collintest(hx)
        % ===== COLINEARITY =====


        % =============== 8: lapse trials
        nc = nc+1;
        if ~isempty(dot_allStart(~sel))

            durations{nc}        = trialDur;
            names{nc}            = 'lapse';
            onsets{nc}           = dot_allStart(~sel);
            orth{nc}             = 0;

            nc = nc+1;

        else
            runList(end) = []
        end



    case 'margRespBlk'    % =====================================================
        %%


        %  ---- resp -----
        resp = r.behav.resp(:, task_runs);
        resp(resp == 1) = -1;
        resp(resp == 3) = 1;

        corResp = r.task.corrResp(:, task_runs);
        corResp(corResp == 1) = -1;
        corResp(corResp == 3) = 1;

        %  ---- RT -----
        rt = r.behav.rt(:, task_runs);
        prevRT = prevRT(:, task_runs);
        nextRT = nextRT(:, task_runs);

        %  ---- Acc -----
        acc = r.behav.acc(:, task_runs);
        prevAcc = prevAcc(:, task_runs);
        nextAcc = nextAcc(:, task_runs);


        %  ---- Dist -----
        dist = r.confs(:, task_runs);
        prevDist = prevDist(:, task_runs);
        nextDist = nextDist(:, task_runs);
        prev2Dist = prev2Dist(:, task_runs);
        next2Dist = next2Dist(:, task_runs);

        %  ---- Targ -----
        targ = r.attendCohs(:, task_runs);
        prevTarg = prevTarg(:, task_runs);
        nextTarg = nextTarg(:, task_runs);
        prev2Targ = prev2Targ(:, task_runs);
        next2Targ = next2Targ(:, task_runs);




        sel = acc==1;


        nc = 1;

        %  =============== 1: response trials

        durations{nc}        = trialDur;
        names{nc}            = 'trial';
        onsets{nc}           = dot_allStart(sel)';
        orth{nc}             = 0;

        pm = 1;



        for bb = 1:size(acc,2)


            blksel = deal(zeros(size(acc)));
            blksel(:, bb) = 1;




            % targ coh
            pmod(nc).name{pm}     = sprintf('targ5R_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)'==5 & corResp(sel)'==1, blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;

            % targ coh
            pmod(nc).name{pm}     = sprintf('targ4R_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)'==4 & corResp(sel)'==1, blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targ2R_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)'==2 & corResp(sel)'==1, blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targ1R_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)'==1 & corResp(sel)'==1, blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;



            % targ coh
            pmod(nc).name{pm}     = sprintf('targ1L_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)'==1 & corResp(sel)'==-1, blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;

            % targ coh
            pmod(nc).name{pm}     = sprintf('targ2L_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)'==2 & corResp(sel)'==-1, blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targ4L_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)'==4 & corResp(sel)'==-1, blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targ5L_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)'==5 & corResp(sel)'==-1, blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;






            % dist coh
            pmod(nc).name{pm}     = sprintf('dist2R_%d', bb);
            pmod(nc).param{pm}    = selz((dist(sel)'==1 & corResp(sel)'==1) | (dist(sel)'==5 & corResp(sel)'==-1), blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('dist1R_%d', bb);
            pmod(nc).param{pm}    = selz((dist(sel)'==2 & corResp(sel)'==1) | (dist(sel)'==4 & corResp(sel)'==-1), blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('dist1L_%d', bb);
            pmod(nc).param{pm}    = selz((dist(sel)'==2 & corResp(sel)'==-1) | (dist(sel)'==4 & corResp(sel)'==1), blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('dist2L_%d', bb);
            pmod(nc).param{pm}    = selz((dist(sel)'==1 & corResp(sel)'==-1) | (dist(sel)'==5 & corResp(sel)'==1), blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


        end




        runList = [1:(length(pmod(1).name)+2)] % keep track of conditions

        % ===== COLINEARITY =====
        xs = []; hx = [];
        for ii = 1:(pm-1)
            xs(:,ii) = pmod(nc).param{ii}';
            hx(:,ii) = conv(pmod(nc).param{ii}', spm_hrf(1.5));
        end

        disp(xs(1:50,1:5))
        statX   = [mean(hx); std(hx)]
        %         corrX   = corr(hx)
        collintest(hx)
        % ===== COLINEARITY =====


        % =============== 8: lapse trials
        nc = nc+1;
        if ~isempty(dot_allStart(~sel))

            durations{nc}        = trialDur;
            names{nc}            = 'lapse';
            onsets{nc}           = dot_allStart(~sel);
            orth{nc}             = 0;

            nc = nc+1;

        else
            runList(end) = []
        end






    case 'FC_PFCl'
        %%


        %  ---- resp -----
        resp = r.behav.resp(:, task_runs);
        resp(resp == 1) = -1;
        resp(resp == 3) = 1;

        corResp = r.task.corrResp(:, task_runs);
        corResp(corResp == 1) = -1;
        corResp(corResp == 3) = 1;

        %  ---- RT -----
        rt = r.behav.rt(:, task_runs);
        prevRT = prevRT(:, task_runs);
        nextRT = nextRT(:, task_runs);

        %  ---- Acc -----
        acc = r.behav.acc(:, task_runs);
        prevAcc = prevAcc(:, task_runs);
        nextAcc = nextAcc(:, task_runs);


        %  ---- Dist -----
        dist = r.confs(:, task_runs);
        prevDist = prevDist(:, task_runs);
        nextDist = nextDist(:, task_runs);
        prev2Dist = prev2Dist(:, task_runs);
        next2Dist = next2Dist(:, task_runs);

        %  ---- Targ -----
        targ = r.attendCohs(:, task_runs);
        prevTarg = prevTarg(:, task_runs);
        nextTarg = nextTarg(:, task_runs);
        prev2Targ = prev2Targ(:, task_runs);
        next2Targ = next2Targ(:, task_runs);




        sel = isfinite(rt);


        nc = 1;

        %  =============== 1: response trials

        durations{nc}        = trialDur;
        names{nc}            = 'trial';
        onsets{nc}           = dot_allStart(sel)';
        orth{nc}             = 0;

        pm = 1;



        for bb = 1:size(acc,2)


            blksel = deal(zeros(size(acc)));
            blksel(:, bb) = 1;


            % targ resp
            pmod(nc).name{pm}     = sprintf('targResp_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist resp
            pmod(nc).name{pm}     = sprintf('distResp_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targCoh_%d', bb);
            pmod(nc).param{pm}    = selz((targ(sel)-3)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('distCoh_%d', bb);
            pmod(nc).param{pm}    = selz(abs(3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('cong_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % RT
            pmod(nc).name{pm}     = sprintf('rt_%d', bb);
            pmod(nc).param{pm}    = selz(rt(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % acc
            pmod(nc).name{pm}     = sprintf('acc_%d', bb);
            pmod(nc).param{pm}    = selz(acc(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;



        end




        runList = [1:(length(pmod(1).name)+2)] % keep track of conditions

        % ===== COLINEARITY =====
        xs = []; hx = [];
        for ii = 1:(pm-1)
            xs(:,ii) = pmod(nc).param{ii}';
            hx(:,ii) = conv(pmod(nc).param{ii}', spm_hrf(1.5));
        end

        disp(xs(1:50,:))
        statX   = [mean(hx); std(hx)]
        corrX   = corr(hx)
        collintest(hx)
        % ===== COLINEARITY =====


        % =============== 8: lapse trials
        nc = nc+1;
        if ~isempty(dot_allStart(~sel))

            durations{nc}        = trialDur;
            names{nc}            = 'lapse';
            onsets{nc}           = dot_allStart(~sel);
            orth{nc}             = 0;

            nc = nc+1;

        else
            runList(end) = []
        end



    case 'FC_IPS-PFCl'
        %%


        %  ---- resp -----
        resp = r.behav.resp(:, task_runs);
        resp(resp == 1) = -1;
        resp(resp == 3) = 1;

        corResp = r.task.corrResp(:, task_runs);
        corResp(corResp == 1) = -1;
        corResp(corResp == 3) = 1;

        %  ---- RT -----
        rt = r.behav.rt(:, task_runs);
        prevRT = prevRT(:, task_runs);
        nextRT = nextRT(:, task_runs);

        %  ---- Acc -----
        acc = r.behav.acc(:, task_runs);
        prevAcc = prevAcc(:, task_runs);
        nextAcc = nextAcc(:, task_runs);


        %  ---- Dist -----
        dist = r.confs(:, task_runs);
        prevDist = prevDist(:, task_runs);
        nextDist = nextDist(:, task_runs);
        prev2Dist = prev2Dist(:, task_runs);
        next2Dist = next2Dist(:, task_runs);

        %  ---- Targ -----
        targ = r.attendCohs(:, task_runs);
        prevTarg = prevTarg(:, task_runs);
        nextTarg = nextTarg(:, task_runs);
        prev2Targ = prev2Targ(:, task_runs);
        next2Targ = next2Targ(:, task_runs);




        sel = isfinite(rt);


        nc = 1;

        %  =============== 1: response trials

        durations{nc}        = trialDur;
        names{nc}            = 'trial';
        onsets{nc}           = dot_allStart(sel)';
        orth{nc}             = 0;

        pm = 1;



        for bb = 1:size(acc,2)


            blksel = deal(zeros(size(acc)));
            blksel(:, bb) = 1;


            % targ resp
            pmod(nc).name{pm}     = sprintf('targResp_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist resp
            pmod(nc).name{pm}     = sprintf('distResp_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targCoh_%d', bb);
            pmod(nc).param{pm}    = selz((targ(sel)-3)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('distCoh_%d', bb);
            pmod(nc).param{pm}    = selz(abs(3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('cong_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % RT
            pmod(nc).name{pm}     = sprintf('rt_%d', bb);
            pmod(nc).param{pm}    = selz(rt(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % acc
            pmod(nc).name{pm}     = sprintf('acc_%d', bb);
            pmod(nc).param{pm}    = selz(acc(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;



        end




        runList = [1:(length(pmod(1).name)+2)] % keep track of conditions

        % ===== COLINEARITY =====
        xs = []; hx = [];
        for ii = 1:(pm-1)
            xs(:,ii) = pmod(nc).param{ii}';
            hx(:,ii) = conv(pmod(nc).param{ii}', spm_hrf(1.5));
        end

        disp(xs(1:50,:))
        statX   = [mean(hx); std(hx)]
        corrX   = corr(hx)
        collintest(hx)
        % ===== COLINEARITY =====


        % =============== 8: lapse trials
        nc = nc+1;
        if ~isempty(dot_allStart(~sel))

            durations{nc}        = trialDur;
            names{nc}            = 'lapse';
            onsets{nc}           = dot_allStart(~sel);
            orth{nc}             = 0;

            nc = nc+1;

        else
            runList(end) = []
        end





        %%


        %  ---- resp -----
        resp = r.behav.resp(:, task_runs);
        resp(resp == 1) = -1;
        resp(resp == 3) = 1;

        corResp = r.task.corrResp(:, task_runs);
        corResp(corResp == 1) = -1;
        corResp(corResp == 3) = 1;

        %  ---- RT -----
        rt = r.behav.rt(:, task_runs);
        prevRT = prevRT(:, task_runs);
        nextRT = nextRT(:, task_runs);

        %  ---- Acc -----
        acc = r.behav.acc(:, task_runs);
        prevAcc = prevAcc(:, task_runs);
        nextAcc = nextAcc(:, task_runs);


        %  ---- Dist -----
        dist = r.confs(:, task_runs);
        prevDist = prevDist(:, task_runs);
        nextDist = nextDist(:, task_runs);
        prev2Dist = prev2Dist(:, task_runs);
        next2Dist = next2Dist(:, task_runs);

        %  ---- Targ -----
        targ = r.attendCohs(:, task_runs);
        prevTarg = prevTarg(:, task_runs);
        nextTarg = nextTarg(:, task_runs);
        prev2Targ = prev2Targ(:, task_runs);
        next2Targ = next2Targ(:, task_runs);




        sel = isfinite(rt);


        nc = 1;

        %  =============== 1: response trials

        durations{nc}        = trialDur;
        names{nc}            = 'trial';
        onsets{nc}           = dot_allStart(sel)';
        orth{nc}             = 0;

        pm = 1;



        for bb = 1:size(acc,2)


            blksel = deal(zeros(size(acc)));
            blksel(:, bb) = 1;


            % targ resp
            pmod(nc).name{pm}     = sprintf('targResp_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist resp
            pmod(nc).name{pm}     = sprintf('distResp_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targCoh_%d', bb);
            pmod(nc).param{pm}    = selz((targ(sel)-3)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('distCoh_%d', bb);
            pmod(nc).param{pm}    = selz(abs(3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('cong_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % RT
            pmod(nc).name{pm}     = sprintf('rt_%d', bb);
            pmod(nc).param{pm}    = selz(rt(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % acc
            pmod(nc).name{pm}     = sprintf('acc_%d', bb);
            pmod(nc).param{pm}    = selz(acc(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;



        end




        runList = [1:(length(pmod(1).name)+2)] % keep track of conditions

        % ===== COLINEARITY =====
        xs = []; hx = [];
        for ii = 1:(pm-1)
            xs(:,ii) = pmod(nc).param{ii}';
            hx(:,ii) = conv(pmod(nc).param{ii}', spm_hrf(1.5));
        end

        disp(xs(1:50,:))
        statX   = [mean(hx); std(hx)]
        corrX   = corr(hx)
        collintest(hx)
        % ===== COLINEARITY =====


        % =============== 8: lapse trials
        nc = nc+1;
        if ~isempty(dot_allStart(~sel))

            durations{nc}        = trialDur;
            names{nc}            = 'lapse';
            onsets{nc}           = dot_allStart(~sel);
            orth{nc}             = 0;

            nc = nc+1;

        else
            runList(end) = []
        end


    case 'FC_IPS'
        %%


        %  ---- resp -----
        resp = r.behav.resp(:, task_runs);
        resp(resp == 1) = -1;
        resp(resp == 3) = 1;

        corResp = r.task.corrResp(:, task_runs);
        corResp(corResp == 1) = -1;
        corResp(corResp == 3) = 1;

        %  ---- RT -----
        rt = r.behav.rt(:, task_runs);
        prevRT = prevRT(:, task_runs);
        nextRT = nextRT(:, task_runs);

        %  ---- Acc -----
        acc = r.behav.acc(:, task_runs);
        prevAcc = prevAcc(:, task_runs);
        nextAcc = nextAcc(:, task_runs);


        %  ---- Dist -----
        dist = r.confs(:, task_runs);
        prevDist = prevDist(:, task_runs);
        nextDist = nextDist(:, task_runs);
        prev2Dist = prev2Dist(:, task_runs);
        next2Dist = next2Dist(:, task_runs);

        %  ---- Targ -----
        targ = r.attendCohs(:, task_runs);
        prevTarg = prevTarg(:, task_runs);
        nextTarg = nextTarg(:, task_runs);
        prev2Targ = prev2Targ(:, task_runs);
        next2Targ = next2Targ(:, task_runs);




        sel = isfinite(rt);


        nc = 1;

        %  =============== 1: response trials

        durations{nc}        = trialDur;
        names{nc}            = 'trial';
        onsets{nc}           = dot_allStart(sel)';
        orth{nc}             = 0;

        pm = 1;



        for bb = 1:size(acc,2)


            blksel = deal(zeros(size(acc)));
            blksel(:, bb) = 1;


            % targ resp
            pmod(nc).name{pm}     = sprintf('targResp_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist resp
            pmod(nc).name{pm}     = sprintf('distResp_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targCoh_%d', bb);
            pmod(nc).param{pm}    = selz((targ(sel)-3)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('distCoh_%d', bb);
            pmod(nc).param{pm}    = selz(abs(3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('cong_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % RT
            pmod(nc).name{pm}     = sprintf('rt_%d', bb);
            pmod(nc).param{pm}    = selz(rt(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % acc
            pmod(nc).name{pm}     = sprintf('acc_%d', bb);
            pmod(nc).param{pm}    = selz(acc(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;



        end




        runList = [1:(length(pmod(1).name)+2)] % keep track of conditions

        % ===== COLINEARITY =====
        xs = []; hx = [];
        for ii = 1:(pm-1)
            xs(:,ii) = pmod(nc).param{ii}';
            hx(:,ii) = conv(pmod(nc).param{ii}', spm_hrf(1.5));
        end

        disp(xs(1:50,:))
        statX   = [mean(hx); std(hx)]
        corrX   = corr(hx)
        collintest(hx)
        % ===== COLINEARITY =====


        % =============== 8: lapse trials
        nc = nc+1;
        if ~isempty(dot_allStart(~sel))

            durations{nc}        = trialDur;
            names{nc}            = 'lapse';
            onsets{nc}           = dot_allStart(~sel);
            orth{nc}             = 0;

            nc = nc+1;

        else
            runList(end) = []
        end


    case 'perfCongBlk'      % =====================================================
        %%


        %  ---- resp -----
        resp = r.behav.resp(:, task_runs);
        resp(resp == 1) = -1;
        resp(resp == 3) = 1;

        corResp = r.task.corrResp(:, task_runs);
        corResp(corResp == 1) = -1;
        corResp(corResp == 3) = 1;

        %  ---- RT -----
        rt = r.behav.rt(:, task_runs);
        prevRT = prevRT(:, task_runs);
        nextRT = nextRT(:, task_runs);

        %  ---- Acc -----
        acc = r.behav.acc(:, task_runs);
        prevAcc = prevAcc(:, task_runs);
        nextAcc = nextAcc(:, task_runs);


        %  ---- Dist -----
        dist = r.confs(:, task_runs);
        prevDist = prevDist(:, task_runs);
        nextDist = nextDist(:, task_runs);
        prev2Dist = prev2Dist(:, task_runs);
        next2Dist = next2Dist(:, task_runs);

        %  ---- Targ -----
        targ = r.attendCohs(:, task_runs);
        prevTarg = prevTarg(:, task_runs);
        nextTarg = nextTarg(:, task_runs);
        prev2Targ = prev2Targ(:, task_runs);
        next2Targ = next2Targ(:, task_runs);




        sel = isfinite(rt);


        nc = 1;

        %  =============== 1: response trials

        durations{nc}        = trialDur;
        names{nc}            = 'trial';
        onsets{nc}           = dot_allStart(sel)';
        orth{nc}             = 0;

        pm = 1;



        for bb = 1:size(acc,2)


            blksel = deal(zeros(size(acc)));
            blksel(:, bb) = 1;


            % targ resp
            pmod(nc).name{pm}     = sprintf('targResp_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist resp
            pmod(nc).name{pm}     = sprintf('distResp_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targCoh_%d', bb);
            pmod(nc).param{pm}    = selz((targ(sel)-3)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('distCoh_%d', bb);
            pmod(nc).param{pm}    = selz(abs(3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('cong_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % RT
            pmod(nc).name{pm}     = sprintf('rt_%d', bb);
            pmod(nc).param{pm}    = selz(rt(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % acc
            pmod(nc).name{pm}     = sprintf('acc_%d', bb);
            pmod(nc).param{pm}    = selz(acc(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;



        end




        runList = [1:(length(pmod(1).name)+2)] % keep track of conditions

        % ===== COLINEARITY =====
        xs = []; hx = [];
        for ii = 1:(pm-1)
            xs(:,ii) = pmod(nc).param{ii}';
            hx(:,ii) = conv(pmod(nc).param{ii}', spm_hrf(1.5));
        end

        disp(xs(1:50,:))
        statX   = [mean(hx); std(hx)]
        corrX   = corr(hx)
        collintest(hx)
        % ===== COLINEARITY =====


        % =============== 8: lapse trials
        nc = nc+1;
        if ~isempty(dot_allStart(~sel))

            durations{nc}        = trialDur;
            names{nc}            = 'lapse';
            onsets{nc}           = dot_allStart(~sel);
            orth{nc}             = 0;

            nc = nc+1;

        else
            runList(end) = []
        end


    case 'feature'    % =====================================================
        %%


        %  ---- resp -----
        resp = r.behav.resp(:, task_runs);
        resp(resp == 1) = -1;
        resp(resp == 3) = 1;

        corResp = r.task.corrResp(:, task_runs);
        corResp(corResp == 1) = -1;
        corResp(corResp == 3) = 1;

        %  ---- RT -----
        rt = r.behav.rt(:, task_runs);
        prevRT = prevRT(:, task_runs);
        nextRT = nextRT(:, task_runs);

        %  ---- Acc -----
        acc = r.behav.acc(:, task_runs);
        prevAcc = prevAcc(:, task_runs);
        nextAcc = nextAcc(:, task_runs);


        %  ---- Dist -----
        dist = r.confs(:, task_runs);
        prevDist = prevDist(:, task_runs);
        nextDist = nextDist(:, task_runs);
        prev2Dist = prev2Dist(:, task_runs);
        next2Dist = next2Dist(:, task_runs);

        %  ---- Targ -----
        targ = r.attendCohs(:, task_runs);
        prevTarg = prevTarg(:, task_runs);
        nextTarg = nextTarg(:, task_runs);
        prev2Targ = prev2Targ(:, task_runs);
        next2Targ = next2Targ(:, task_runs);




        sel = isfinite(rt);


        nc = 1;

        %  =============== 1: response trials

        durations{nc}        = trialDur;
        names{nc}            = 'trial';
        onsets{nc}           = dot_allStart(sel)';
        orth{nc}             = 0;

        pm = 1;


        % targ resp
        pmod(nc).name{pm}     = sprintf('targResp');
        pmod(nc).param{pm}    = nanz(targ(sel)' .* corResp(sel)');
        pmod(nc).poly{pm}     = 1;
        pm = pm+1;


        % dist resp
        pmod(nc).name{pm}     = sprintf('distResp');
        pmod(nc).param{pm}    = nanz((3-dist(sel))' .* corResp(sel)');
        pmod(nc).poly{pm}     = 1;
        pm = pm+1;



        % targ coh
        pmod(nc).name{pm}     = sprintf('targCoh');
        pmod(nc).param{pm}    = nanz((targ(sel)-3)');
        pmod(nc).poly{pm}     = 1;
        pm = pm+1;


        % dist coh
        pmod(nc).name{pm}     = sprintf('distCoh');
        pmod(nc).param{pm}    = nanz(abs(3-dist(sel))');
        pmod(nc).poly{pm}     = 1;
        pm = pm+1;


        % dist coh
        pmod(nc).name{pm}     = sprintf('cong');
        pmod(nc).param{pm}    = nanz((3-dist(sel))');
        pmod(nc).poly{pm}     = 1;
        pm = pm+1;





        runList = [1:(length(pmod(1).name)+2)] % keep track of conditions

        % ===== COLINEARITY =====
        xs = []; hx = [];
        for ii = 1:(pm-1)
            xs(:,ii) = pmod(nc).param{ii}';
            hx(:,ii) = conv(pmod(nc).param{ii}', spm_hrf(1.5));
        end

        disp(xs(1:50,1:5))
        statX   = [mean(hx); std(hx)]
        %         corrX   = corr(hx)
        collintest(hx)
        % ===== COLINEARITY =====


        % =============== 8: lapse trials
        nc = nc+1;
        if ~isempty(dot_allStart(~sel))

            durations{nc}        = trialDur;
            names{nc}            = 'lapse';
            onsets{nc}           = dot_allStart(~sel);
            orth{nc}             = 0;

            nc = nc+1;

        else
            runList(end) = []
        end



    case 'perfBlk-comb'    % =====================================================
        %%


        %  ---- resp -----
        resp = r.behav.resp(:, task_runs);
        resp(resp == 1) = -1;
        resp(resp == 3) = 1;

        corResp = r.task.corrResp(:, task_runs);
        corResp(corResp == 1) = -1;
        corResp(corResp == 3) = 1;

        %  ---- RT -----
        rt = r.behav.rt(:, task_runs);
        prevRT = prevRT(:, task_runs);
        nextRT = nextRT(:, task_runs);

        %  ---- Acc -----
        acc = r.behav.acc(:, task_runs);
        prevAcc = prevAcc(:, task_runs);
        nextAcc = nextAcc(:, task_runs);


        %  ---- Dist -----
        dist = r.confs(:, task_runs);
        prevDist = prevDist(:, task_runs);
        nextDist = nextDist(:, task_runs);
        prev2Dist = prev2Dist(:, task_runs);
        next2Dist = next2Dist(:, task_runs);

        %  ---- Targ -----
        targ = r.attendCohs(:, task_runs);
        prevTarg = prevTarg(:, task_runs);
        nextTarg = nextTarg(:, task_runs);
        prev2Targ = prev2Targ(:, task_runs);
        next2Targ = next2Targ(:, task_runs);




        sel = (acc==1) & (cumsum(acc) <= 45);


        nc = 1;

        %  =============== 1: response trials

        durations{nc}        = trialDur;
        names{nc}            = 'trial';
        onsets{nc}           = dot_allStart(sel)';
        orth{nc}             = 0;

        pm = 1;



        for bb = 1:size(acc,2)


            blksel = deal(zeros(size(acc)));
            blksel(:, bb) = 1;


            % targ resp
            pmod(nc).name{pm}     = sprintf('targResp_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist resp
            pmod(nc).name{pm}     = sprintf('distResp_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targCoh_%d', bb);
            pmod(nc).param{pm}    = selz((targ(sel)-3)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('distCoh_%d', bb);
            pmod(nc).param{pm}    = selz(abs(3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('cong_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;



        end


        % dist coh
        pmod(nc).name{pm}     = sprintf('rt');
        pmod(nc).param{pm}    = nanz(rt(sel)');
        pmod(nc).poly{pm}     = 1;
        pm = pm+1;





        runList = [1:(length(pmod(1).name)+2)] % keep track of conditions

        % ===== COLINEARITY =====
        xs = []; hx = [];
        for ii = 1:(pm-1)
            xs(:,ii) = pmod(nc).param{ii}';
            hx(:,ii) = conv(pmod(nc).param{ii}', spm_hrf(1.5));
        end

        disp(xs(1:50,1:5))
        statX   = [mean(hx); std(hx)]
        %         corrX   = corr(hx)
        collintest(hx)
        % ===== COLINEARITY =====


        % =============== 8: lapse trials
        nc = nc+1;
        if ~isempty(dot_allStart(~sel))

            durations{nc}        = trialDur;
            names{nc}            = 'lapse';
            onsets{nc}           = dot_allStart(~sel);
            orth{nc}             = 0;

            nc = nc+1;

        else
            runList(end) = []
        end



    case 'featureBlk-comb'    % =====================================================
        %%


        %  ---- resp -----
        resp = r.behav.resp(:, task_runs);
        resp(resp == 1) = -1;
        resp(resp == 3) = 1;

        corResp = r.task.corrResp(:, task_runs);
        corResp(corResp == 1) = -1;
        corResp(corResp == 3) = 1;

        %  ---- RT -----
        rt = r.behav.rt(:, task_runs);
        prevRT = prevRT(:, task_runs);
        nextRT = nextRT(:, task_runs);

        %  ---- Acc -----
        acc = r.behav.acc(:, task_runs);
        prevAcc = prevAcc(:, task_runs);
        nextAcc = nextAcc(:, task_runs);


        %  ---- Dist -----
        dist = r.confs(:, task_runs);
        prevDist = prevDist(:, task_runs);
        nextDist = nextDist(:, task_runs);
        prev2Dist = prev2Dist(:, task_runs);
        next2Dist = next2Dist(:, task_runs);

        %  ---- Targ -----
        targ = r.attendCohs(:, task_runs);
        prevTarg = prevTarg(:, task_runs);
        nextTarg = nextTarg(:, task_runs);
        prev2Targ = prev2Targ(:, task_runs);
        next2Targ = next2Targ(:, task_runs);




        sel = acc==1 & cumsum(acc)<45;


        nc = 1;

        %  =============== 1: response trials

        durations{nc}        = trialDur;
        names{nc}            = 'trial';
        onsets{nc}           = dot_allStart(sel)';
        orth{nc}             = 0;

        pm = 1;



        for bb = 1:size(acc,2)


            blksel = deal(zeros(size(acc)));
            blksel(:, bb) = 1;


            % targ resp
            pmod(nc).name{pm}     = sprintf('targResp_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist resp
            pmod(nc).name{pm}     = sprintf('distResp_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targCoh_%d', bb);
            pmod(nc).param{pm}    = selz((targ(sel)-3)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('distCoh_%d', bb);
            pmod(nc).param{pm}    = selz(abs(3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('cong_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;



        end




        runList = [1:(length(pmod(1).name)+2)] % keep track of conditions

        % ===== COLINEARITY =====
        xs = []; hx = [];
        for ii = 1:(pm-1)
            xs(:,ii) = pmod(nc).param{ii}';
            hx(:,ii) = conv(pmod(nc).param{ii}', spm_hrf(1.5));
        end

        disp(xs(1:50,1:5))
        statX   = [mean(hx); std(hx)]
        %         corrX   = corr(hx)
        collintest(hx)
        % ===== COLINEARITY =====


        % =============== 8: lapse trials
        nc = nc+1;
        if ~isempty(dot_allStart(~sel))

            durations{nc}        = trialDur;
            names{nc}            = 'lapse';
            onsets{nc}           = dot_allStart(~sel);
            orth{nc}             = 0;

            nc = nc+1;

        else
            runList(end) = []
        end




    case 'perfBlk'    % =====================================================
        %%


        %  ---- resp -----
        resp = r.behav.resp(:, task_runs);
        resp(resp == 1) = -1;
        resp(resp == 3) = 1;

        corResp = r.task.corrResp(:, task_runs);
        corResp(corResp == 1) = -1;
        corResp(corResp == 3) = 1;

        %  ---- RT -----
        rt = r.behav.rt(:, task_runs);
        prevRT = prevRT(:, task_runs);
        nextRT = nextRT(:, task_runs);

        %  ---- Acc -----
        acc = r.behav.acc(:, task_runs);
        prevAcc = prevAcc(:, task_runs);
        nextAcc = nextAcc(:, task_runs);


        %  ---- Dist -----
        dist = r.confs(:, task_runs);
        prevDist = prevDist(:, task_runs);
        nextDist = nextDist(:, task_runs);
        prev2Dist = prev2Dist(:, task_runs);
        next2Dist = next2Dist(:, task_runs);

        %  ---- Targ -----
        targ = r.attendCohs(:, task_runs);
        prevTarg = prevTarg(:, task_runs);
        nextTarg = nextTarg(:, task_runs);
        prev2Targ = prev2Targ(:, task_runs);
        next2Targ = next2Targ(:, task_runs);




        sel = isfinite(rt);


        nc = 1;

        %  =============== 1: response trials

        durations{nc}        = trialDur;
        names{nc}            = 'trial';
        onsets{nc}           = dot_allStart(sel)';
        orth{nc}             = 0;

        pm = 1;



        for bb = 1:size(acc,2)


            blksel = deal(zeros(size(acc)));
            blksel(:, bb) = 1;


            % targ resp
            pmod(nc).name{pm}     = sprintf('targResp_%d', bb);
            pmod(nc).param{pm}    = selz(targ(sel)' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist resp
            pmod(nc).name{pm}     = sprintf('distResp_%d', bb);
            pmod(nc).param{pm}    = selz((3-dist(sel))' .* corResp(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % targ coh
            pmod(nc).name{pm}     = sprintf('targCoh_%d', bb);
            pmod(nc).param{pm}    = selz((targ(sel)-3)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % dist coh
            pmod(nc).name{pm}     = sprintf('distCoh_%d', bb);
            pmod(nc).param{pm}    = selz(abs(3-dist(sel))', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % RT
            pmod(nc).name{pm}     = sprintf('rt_%d', bb);
            pmod(nc).param{pm}    = selz(rt(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;


            % acc
            pmod(nc).name{pm}     = sprintf('acc_%d', bb);
            pmod(nc).param{pm}    = selz(acc(sel)', blksel(sel)');
            pmod(nc).poly{pm}     = 1;
            pm = pm+1;



            %              % RT
            %             pmod(nc).name{pm}     = sprintf('targRt_%d', bb);
            %             pmod(nc).param{pm}    = selz(rt(sel)', blksel(sel)') .* selz((targ(sel)-3)', blksel(sel)');
            %             pmod(nc).poly{pm}     = 1;
            %             pm = pm+1;
            %
            %
            %              % RT
            %             pmod(nc).name{pm}     = sprintf('targAcc_%d', bb);
            %             pmod(nc).param{pm}    = selz(acc(sel)', blksel(sel)') .* selz((targ(sel)-3)', blksel(sel)');
            %             pmod(nc).poly{pm}     = 1;
            %             pm = pm+1;
            %
            %
            %             % acc
            %             pmod(nc).name{pm}     = sprintf('distRt_%d', bb);
            %             pmod(nc).param{pm}    = selz(rt(sel)', blksel(sel)') .* selz((3-dist(sel))', blksel(sel)');
            %             pmod(nc).poly{pm}     = 1;
            %             pm = pm+1;
            %
            %
            %              % acc
            %             pmod(nc).name{pm}     = sprintf('distAcc_%d', bb);
            %             pmod(nc).param{pm}    = selz(acc(sel)', blksel(sel)') .* selz((3-dist(sel))', blksel(sel)');
            %             pmod(nc).poly{pm}     = 1;
            %             pm = pm+1;




        end




        runList = [1:(length(pmod(1).name)+2)] % keep track of conditions

        % ===== COLINEARITY =====
        xs = []; hx = [];
        for ii = 1:(pm-1)
            xs(:,ii) = pmod(nc).param{ii}';
            hx(:,ii) = conv(pmod(nc).param{ii}', spm_hrf(1.5));
        end

        disp(xs(1:50,:))
        statX   = [mean(hx); std(hx)]
        corrX   = corr(hx)
        collintest(hx)
        % ===== COLINEARITY =====


        % =============== 8: lapse trials
        nc = nc+1;
        if ~isempty(dot_allStart(~sel))

            durations{nc}        = trialDur;
            names{nc}            = 'lapse';
            onsets{nc}           = dot_allStart(~sel);
            orth{nc}             = 0;

            nc = nc+1;

        else
            runList(end) = []
        end


    
    
    otherwise 'NO NAME'


end




% save mat
save(sprintf('%s/%s', cond_dir, analysis),...
    'names', 'onsets', 'durations', 'pmod', 'orth')



%% === make confound regressors

R = [];

csfmx = [];
runmx = [];
hpfmx = [];
all_fd = [];
for rr = 1:nRuns

    % load data
    cs = tdfread(fullfile(cl(rr).folder, cl(rr).name));

    % framewise displacement
    fd = [0; str2num(cs.framewise_displacement(2:end,:))];
    all_fd = [all_fd; fd];



    % HPF
    if any(regName{1}) % if MOT

        % cosine HPF
        cos0 = cs.cosine00;

        hpfmx = blkdiag(hpfmx, [cos0]);

    elseif any(regName{2}) % if COMBINED

        if any(task_runs(rr) == colSel)

            % cosine HPF
            cos0 = cs.cosine00;
            cos1 = cs.cosine01;
            cos2 = cs.cosine02;
            cos3 = cs.cosine03;
            cos4 = cs.cosine04;
            cos5 = cs.cosine05;

            hpfmx = blkdiag(hpfmx, [cos0, cos1, cos2, cos3, cos4, cos5]);

        else

            % cosine HPF
            cos0 = cs.cosine00;

            hpfmx = blkdiag(hpfmx, [cos0]);

        end

    elseif any(regName{3}) % if localizer

        % cosine HPF
        cos0 = cs.cosine00;
        cos1 = cs.cosine01;

        hpfmx = blkdiag(hpfmx, [cos0, cos1]);


    else

        % cosine HPF
        cos0 = cs.cosine00;
        cos1 = cs.cosine01;
        cos2 = cs.cosine02;
        cos3 = cs.cosine03;
        cos4 = cs.cosine04;
        cos5 = cs.cosine05;

        hpfmx = blkdiag(hpfmx, [cos0, cos1, cos2, cos3, cos4, cos5]);

    end



    % run intercepts
    runmx = blkdiag(runmx, ones(size(cs.trans_x)));



    %     % CSF-WM
    csfwm   = center(cs.csf_wm);
    % confound regressors
    R = [R; ...
        csfwm, ...
        ];

    %     csfmx = blkdiag(csfmx, zscore(cs.csf_wm));
    %     names{rr} = sprintf('csfwm_%d', rr);


end


%  R = [R, csfmx];


names = { ...
    'csfwm', ...
    };


switch cvi

    case 'fast'

        names{end+1} = 'fd';
        R = [R, all_fd];


    case 'wls'

        % HPF (125s)
        for rr = 1:nRuns

            if any(regName{1}) % if MOT

                names{end+1} = sprintf('hpf0_%d',rr);


            elseif any(regName{2}) % if COMBINED

                if any(task_runs(rr) == colSel)

                    names{end+1} = sprintf('hpf0_%d',rr);
                    names{end+1} = sprintf('hpf1_%d',rr);
                    names{end+1} = sprintf('hpf2_%d',rr);
                    names{end+1} = sprintf('hpf3_%d',rr);
                    names{end+1} = sprintf('hpf4_%d',rr);
                    names{end+1} = sprintf('hpf5_%d',rr);

                else
                    names{end+1} = sprintf('hpf0_%d',rr);

                end

            elseif any(regName{3}) % if localizer

                names{end+1} = sprintf('hpf0_%d',rr);
                names{end+1} = sprintf('hpf1_%d',rr);


            else

                names{end+1} = sprintf('hpf0_%d',rr);
                names{end+1} = sprintf('hpf1_%d',rr);
                names{end+1} = sprintf('hpf2_%d',rr);
                names{end+1} = sprintf('hpf3_%d',rr);
                names{end+1} = sprintf('hpf4_%d',rr);
                names{end+1} = sprintf('hpf5_%d',rr);

            end

        end
        R = [R, hpfmx];

end


% n-1 run intercepts
for rr = 2:nRuns

    names{end+1} = sprintf('run_%d',rr);

end
R = [R, runmx(:,2:end)];




switch name

    case 'FC_PFCl'
        %%
        fc_dir =  sprintf('%s/spm-data/FC/data', root_dir)

        clear D
        D = load(fullfile(fc_dir, sprintf('%d_perfCongBlk_PFCl', ptNum)));


        % prepend FC
        ppiR = [];
        rc = 1;
        for rr = 1:nRuns
            ppiR = [ppiR, selz(D.pc1, runmx(:,rr))];
            ppiName{rc} = sprintf('PFCl_%d', rr); rc = rc+1;
        end
        R = [ppiR, R];
        names = [ppiName, names];

        % save
        save(sprintf('%s/confound', cfd_dir),...
            'R', 'names');

        % update list
        runList = [runList, 1:length(ppiName), zeros(size(names))];
        condList = [condList, runList];
        all_R = [all_R; R];



    case 'FC_IPS-PFCl'
        %%
        fc_dir =  sprintf('%s/spm-data/FC/data', root_dir)

        clear D
        D{1} = load(fullfile(fc_dir, sprintf('%d_perfCongBlk_IPS', ptNum)));
        D{2} = load(fullfile(fc_dir, sprintf('%d_perfCongBlk_PFCl', ptNum)));

        reg_corr = corr(D{1}.pc1, D{2}.pc1)

        % prepend FC
        ppiR = [];
        rc = 1;
        for rr = 1:nRuns
            ppiR = [ppiR, selz(D{1}.pc1, runmx(:,rr))];
            ppiName{rc} = sprintf('IPS_%d', rr); rc = rc+1;
        end
        for rr = 1:nRuns
            ppiR = [ppiR, selz(D{2}.pc1, runmx(:,rr))];
            ppiName{rc} = sprintf('PFCl_%d', rr); rc = rc+1;
        end
        R = [ppiR, R];
        names = [ppiName, names];

        % save
        save(sprintf('%s/confound', cfd_dir),...
            'R', 'names', 'reg_corr');

        % update list
        runList = [runList, 1:length(ppiName), zeros(size(names))];
        condList = [condList, runList];
        all_R = [all_R; R];



    case 'FC_IPS'
        %%
        fc_dir =  sprintf('%s/spm-data/FC/data', root_dir)

        clear D
        D = load(fullfile(fc_dir, sprintf('%d_perfCongBlk_IPS', ptNum)));


        % prepend FC
        ppiR = [];
        rc = 1;
        for rr = 1:nRuns
            ppiR = [ppiR, selz(D.pc1, runmx(:,rr))];
            ppiName{rc} = sprintf('IPS_%d', rr); rc = rc+1;
        end
        R = [ppiR, R];
        names = [ppiName, names];

        % save
        save(sprintf('%s/confound', cfd_dir),...
            'R', 'names');

        % update list
        runList = [runList, 1:length(ppiName), zeros(size(names))];
        condList = [condList, runList];
        all_R = [all_R; R];





    otherwise

        % save & update list
        save(sprintf('%s/confound', cfd_dir),...
            'R', 'names');


        runList = [runList, zeros(size(names))]
        condList = [condList, runList];

        all_R = [all_R; R];

end



save(sprintf('%s/%s_condList', cond_dir, analysis), 'condList', 'task_runs', 'taskMot', 'taskCol');



%% plot motion
mkdir(sprintf('%s/spm-data/motion', root_dir));

f=figure; hold on;
f.Position = [100 100 1300 1000];

% Rs = 150:150:900;


% % % translation % % %
subplot(3,1,1); hold on;
plot(all_R(:,ismember(names, {'tx', 'ty', 'tz'})), '-', 'LineWidth', 2);
% for rr = 1:length(Rs)
%     xline(Rs(rr), '-k', 'LineWidth', .75);
% end
legend({'tx', 'ty', 'tz'}, 'Location', 'northeastoutside')
title(sprintf('translation - sub-%d', ptNum))
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlabel('TR')


% % % rotation % % %
subplot(3,1,2); hold on;
plot(all_R(:,ismember(names, {'rx', 'ry', 'rz'})), '-', 'LineWidth', 2);
% for rr = 1:length(Rs)
%     xline(Rs(rr), '-k', 'LineWidth', .75);
% end
legend({'rx', 'ry', 'rz'}, 'Location', 'northeastoutside')
title(sprintf('rotation - sub-%d', ptNum))
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlabel('TR')


% % % displacement % % %
subplot(3,1,3); hold on;
plot(all_fd, '-', 'LineWidth', 2);
% for rr = 1:length(Rs)
%     xline(Rs(rr), '-k', 'LineWidth', .75);
% end
legend({'fd'}, 'Location', 'northeastoutside')
title(sprintf('displacement - sub-%d', ptNum))
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlabel('TR')


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




% % % save & close % % %
saveas(f,sprintf('%s/spm-data/motion/%s_sub-%d.png', root_dir, movtName, ptNum))

close(f);


%% === get masks

ml = dir(sprintf('%s/anat/%s', pt_dir, mask_filt));
mask_file = fullfile(ml(end).folder, ml(end).name)

ml = dir(sprintf('%s/anat/%s', pt_dir, grey_filt));
grey_file = fullfile(ml(end).folder, ml(end).name)



%% === BATCH ====
%  ==============

spm_jobman('initcfg')

% specify model for each run
for rr = 1:nRuns

    % get frames
    runFile = fullfile(fl(rr).folder, fl(rr).name)
    matlabbatch{rr}.spm.util.exp_frames.files = {runFile};
    matlabbatch{rr}.spm.util.exp_frames.frames = Inf;


    % add scan/cond info to session
    matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess.scans(rr) = ...
        cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{rr}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));

end


matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});




switch name

    case 'cDistPPI' % no conds if PPI

        matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess.multi = {''};

    otherwise

        matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess.multi = ...
            {sprintf('%s/%s.mat', cond_dir, analysis)};
end


matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess.multi_reg = {sprintf('%s/confound.mat', cfd_dir)};
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.sess.hpf = inf; % remove HPF, include cosine basis for each block





% specify general model features
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.dir = {save_dir};
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.timing.units = 'secs';
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.timing.RT = data.param.TR;
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.timing.fmri_t = 16;
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.timing.fmri_t0 = 8;

matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.fact = struct('name', {}, 'levels', {});
switch name
    case 'perfCongBlkD'
        matlabbatch{nRuns+1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
    otherwise
        matlabbatch{nRuns+1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
end
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.volt = 1;
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.global = 'None';
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.mthresh = 0.8;
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.mask = {mask_file};
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.cvi = cvi;
matlabbatch{nRuns+1}.spm.tools.rwls.fmri_rwls_spec.cvi_mask = {grey_file};



% estimate model
matlabbatch{nRuns+2}.spm.tools.rwls.fmri_rwls_est.spmmat(1) = ...
    cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{nRuns+1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{nRuns+2}.spm.tools.rwls.fmri_rwls_est.method.Classical = 1;
% matlabbatch{nRuns+2}.spm.tools.rwls.fmri_rwls_est.write_residuals = 1;




%% === RUN

tic
spm_jobman('run', matlabbatch);
toc


%% === Write Residuals
if strcmp(confound, 'resid')

    disp('calculate residual')

    tic
    load(fullfile(save_dir, 'SPM.mat'));
    spm_write_residuals(SPM, NaN);
    toc

end


end