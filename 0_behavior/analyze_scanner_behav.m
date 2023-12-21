%% analyze pilot data
clear; clc;

pts = 8004:8032

npt = length(pts)



%% load data

vec = @(x) x(:);
center = @(x) x-nanmean(x);

[colTarg, colTarg1, colDist, colRT, colCRT, colAcc, colPt,...
    colAcc1, motAcc1,...
    colDist1, motDist1,...
    colDistCoh,...
    motTarg, motDist, motRT, motAcc, motPt, ...
    colRR, colStim, colAltStim, colPES,...
    colITI, motITI,...
    colResp, colDur, colLapse,...
    colTrial, motTrial,...
    colRuns] = deal([]);


clear all_acc all_lapse demo_gender demo_genderBin demo_age;

for pp = 1:npt

    dr = dir(sprintf('./*dots%d*',pts(pp)));

    d = load(fullfile(dr.folder, dr.name));
    s = d.data.results.rand.session;

    colBlk = sum(isfinite(s.behav.rt)) > 50;

    switch pts(pp)
        case 8010
            colBlk(10) = 0;
        case 8019
            colBlk(10) = 0;
    end

    GoodBlk(pp) = sum(colBlk)./6;

    if pts(pp) == 8006

        demo_gender{pp} = 'M';
        demo_genderBin(pp) = 0;
        demo_age(pp) = 20;

    else

        demo_gender{pp} = d.data.ptInfo{3};
        demo_genderBin(pp) = strcmp(d.data.ptInfo(3), 'F');
        demo_age(pp) = d.data.ptInfo{2};

    end


    % trial
    colTrial = [colTrial; vec(cumsum(isfinite(s.behav.acc(:, colBlk))))];
    motTrial = [motTrial; vec(cumsum(isfinite(s.behav.acc(:,~colBlk))))];


    %     % block
    %     colTrial = [colTrial; vec(cumsum(isfinite(s.behav.acc(:, colBlk))))];
    %     motTrial = [motTrial; vec(cumsum(isfinite(s.behav.acc(:,~colBlk))))];


    % ITI
    colITI = [colITI; vec([nan(1, sum(colBlk)); (s.task.ITI_len(1:end-1, colBlk))])];
    motITI = [motITI; vec([nan(1, sum(~colBlk)); (s.task.ITI_len(1:end-1, ~colBlk))])];

    % targ
    colTarg = [colTarg; vec(s.attendCohs(:, colBlk)-3)];
    motTarg = [motTarg; vec(s.attendCohs(:,~colBlk)-3)];

    colTarg1 = [colTarg1; vec([nan(1, sum(colBlk)); s.attendCohs(1:end-1,  colBlk)-3])];


    % dist
    colDist = [colDist; vec(3-s.confs(:,  colBlk))];
    colDistCoh = [colDistCoh; center(abs(vec(3-s.confs(:,  colBlk))))];
    motDist = [motDist; -center(discretize(vec((s.confs(:, ~colBlk))),5))];

    % prev dist
    colDist1 = [colDist1; vec([nan(1, sum(colBlk)); 3-s.confs(1:end-1,  colBlk)])];
    motDist1 = [motDist1; vec([nan(1, sum(~colBlk)); -center(s.confs(1:end-1, ~colBlk))])];

    % acc
    colAcc = [colAcc; vec(s.behav.acc(:, colBlk))];
    motAcc = [motAcc; vec(s.behav.acc(:, ~colBlk))];

    colAcc1 = [colAcc1; vec([nan(1, sum(colBlk)); s.behav.acc(1:end-1, colBlk)])];
    motAcc1 = [motAcc1; vec([nan(1, sum(~colBlk)); s.behav.acc(1:end-1, ~colBlk)])];


    colPES = [colPES; vec([s.behav.acc(2:end, colBlk); nan(1, sum(colBlk))] - [nan(1, sum(colBlk)); s.behav.acc(1:end-1, colBlk)])];

    % lapse
    colLapse = [colLapse; vec(isfinite(s.behav.rt(:, colBlk)))];



    % resp repeat
    colRR = [colRR; center(vec([nan(1, sum(colBlk)); abs(diff(s.task.corrResp(:, colBlk)))==0]))];

    % rt
    colRT = [colRT; vec(s.behav.rt(:, colBlk))];
    motRT = [motRT; vec(s.behav.rt(:, ~colBlk))];

    colCRT = [colCRT; center(vec(s.behav.rt(:, colBlk)))];



    % color
    colStim = [colStim; vec(s.task.color(:, colBlk))];
    colAltStim = [colAltStim; vec(s.task.altCol(:, colBlk))];

    % resp
    colResp = [colResp; vec(s.behav.resp(:, colBlk)==3)];


    % dur
    colDur  = [colDur; vec(s.time.trOffset(:, colBlk) - s.time.trOnset(:, colBlk))];

    % block
    colRuns = [colRuns; vec(isfinite(s.behav.acc(:, colBlk)) .* [1:sum(colBlk)])];


    % pt
    colPt = [colPt; vec(pts(pp) + 0*s.behav.acc(:, colBlk))];
    motPt = [motPt; vec(pts(pp) + 0*s.behav.acc(:, ~colBlk))];


    meanResp = mean(isfinite(vec(s.behav.rt(:, colBlk))))
    meanAcc = mean(vec(s.behav.acc(:, colBlk)))


    all_acc(pp) = meanAcc; 
    all_lapse(pp) = 1-meanResp;




    %     rval = corr(vec(s.time.trOffset(:, colBlk) - s.time.trOnset(:, colBlk)), vec(s.behav.rt(:, colBlk)), 'rows', 'pairwise')
    %
    %             figure; hold on;
    %             histogram(vec(s.time.trOffset - s.time.trOnset))
    %             title(pts(pp));

end

accAverage = [mean(all_acc), min(all_acc), max(all_acc)]
lapseAverage = [mean(all_lapse), min(all_lapse), max(all_lapse)]
female = mean(demo_genderBin)*100
age = [mean(demo_age), std(demo_age)]

tcol = table;
tcol.targ = colTarg;
tcol.targ1 = colTarg1;
tcol.dist = colDist;
tcol.distCoh = colDistCoh;
tcol.dist1 = colDist1;
tcol.sqdist = abs(colDist);
tcol.sqdist1 = colDist1.^2;
tcol.miss = isnan(colRT);
tcol.acc = colAcc;
tcol.acc1 = colAcc1;
tcol.rt = colRT;
tcol.crt = colCRT;
tcol.resp = colResp;
tcol.lapse = colLapse;
tcol.dur = colDur;
tcol.lrt = log(colRT);
tcol.pes = colPES;
tcol.RR = colRR;
tcol.ITI = colITI;
tcol.trial = colTrial;
tcol.stim = colStim;
tcol.altStim = colAltStim;
tcol.pt = colPt;
tcol.run = colRuns;
tcol.crun = colRuns-3;



tmot = table;
tmot.targ = motTarg;
tmot.dist = motDist;
tmot.dist1 = motDist1;
tmot.miss = isnan(motRT);
tmot.acc = motAcc;
tmot.rt = motRT;
tmot.lrt = log(motRT);
tmot.ITI = motITI;
tmot.trial = motTrial;
tmot.pt = motPt;



% set up regression
rtbl = tcol;
rtbl.ITI = categorical(rtbl.ITI, 'Ordinal', 1);
rtbl.resp = categorical(rtbl.resp);
selCol = isfinite(rtbl.pt) & rtbl.acc1 == 1;


rtblMot = tmot;
rtblMot.targ = center(tmot.targ);
selMot = isfinite(rtbl.pt) & rtbl.acc1 == 1;




%% plot AC vs AM

figure; hold on


[f,x] = ksdensity(tcol.rt(tcol.acc == 1 & tcol.rt>.200));
plot(x,f, '-b', 'LineWidth', 2)

[f,x] = ksdensity(tmot.rt(tmot.acc == 1  & tmot.rt>.200));
plot(x,f, '-r', 'LineWidth', 2)

set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlabel('Reaction Time (s)')
ylabel('PDF')



median(tcol.rt(tcol.acc == 1))
median(tmot.rt(tmot.acc == 1))



%% ========== run regression ============



%% within-trial

selCol = isfinite(rtbl.pt);


% no interaction
within_rt = fitlme(rtbl, 'lrt ~ (targ+dist+distCoh)  + (1+(targ+dist+distCoh)  | pt)', 'Exclude', ~(selCol & rtbl.acc == 1),...
    'FitMethod', 'REML', 'DummyVarCoding', 'effects');

[~,~,ffx_rt] = fixedEffects(within_rt, 'DFMethod', 'Satterthwaite')
[~,~,rfx_rt] = randomEffects(within_rt);
ffx_rt.D = ffx_rt.tStat ./ sqrt(ffx_rt.DF+1)
ffx_rt.Lower = round(ffx_rt.Lower,3, 'significant')
ffx_rt.Upper = round(ffx_rt.Upper,3, 'significant')


within_acc = fitglme(rtbl, 'acc ~ (targ+dist+distCoh) + (1+(targ+dist+distCoh) | pt)', 'Distribution', 'binomial', 'Exclude', ~(selCol & isfinite(rtbl.rt)),...
    'FitMethod', 'Laplace', 'DummyVarCoding', 'effects');

[~,~,ffx_acc] = fixedEffects(within_acc);
ffx_acc.DF = repmat(npt-1, [3,1]);
ffx_acc.pValue = tcdf(abs(ffx_acc.tStat), npt-1, 'upper')*2
ffx_acc.D = ffx_acc.tStat ./ sqrt(npt)
ffx_acc.Lower = round(ffx_acc.Lower,3, 'significant')
ffx_acc.Upper = round(ffx_acc.Upper,3, 'significant')

[~,~,rfx_acc] = randomEffects(within_acc);




[targ_r, targ_p] = corr(rfx_rt.Estimate(ismember(rfx_rt.Name, 'targ')), rfx_acc.Estimate(ismember(rfx_acc.Name, 'targ')))
[dist_r, dist_p] = corr(rfx_rt.Estimate(ismember(rfx_rt.Name, 'dist')), rfx_acc.Estimate(ismember(rfx_acc.Name, 'dist')))


baseRT = rfx_rt.Estimate(ismember(rfx_rt.Name, '(Intercept)'));
targRT = rfx_rt.Estimate(ismember(rfx_rt.Name, 'targ'));
distRT = rfx_rt.Estimate(ismember(rfx_rt.Name, 'dist'));
baseAcc = rfx_acc.Estimate(ismember(rfx_acc.Name, '(Intercept)'));
targAcc = rfx_acc.Estimate(ismember(rfx_acc.Name, 'targ'));
distAcc = rfx_acc.Estimate(ismember(rfx_acc.Name, 'dist'));


save('targdist_rfx', 'baseRT', 'targRT', 'distRT', 'baseAcc', 'targAcc', 'distAcc');



% with interaction
within_rtX = fitlme(rtbl, 'lrt ~ (targ*dist)  + (1+(targ*dist)  | pt)', 'Exclude', ~(selCol & rtbl.acc == 1),...
    'FitMethod', 'REML', 'DummyVarCoding', 'effects');
[~,~,ffx_rtX] = fixedEffects(within_rtX, 'DFMethod', 'Satterthwaite')
ffx_rtX.D = ffx_rtX.tStat ./ sqrt(npt)
ffx_rtX.Lower = round(ffx_rtX.Lower,3, 'significant')
ffx_rtX.Upper = round(ffx_rtX.Upper,3, 'significant')

within_accX = fitglme(rtbl, 'acc ~ (targ*dist) + (1+(targ*dist) | pt)', 'Distribution', 'binomial', 'Exclude', ~(selCol & isfinite(rtbl.rt)),...
    'FitMethod', 'Laplace', 'DummyVarCoding', 'effects');


[~,~,ffx_accX] = fixedEffects(within_accX);
ffx_accX.DF = repmat(npt-1, [4,1]);
ffx_accX.pValue = tcdf(abs(ffx_accX.tStat), npt-1, 'upper')*2

ffx_accX.D = ffx_accX.tStat ./ sqrt(npt)
ffx_accX.Lower = round(ffx_accX.Lower,3, 'significant')
ffx_accX.Upper = round(ffx_accX.Upper,3, 'significant')


% model comparision 
within_rt.ModelCriterion.AIC - within_rtX.ModelCriterion.AIC
within_acc.ModelCriterion.AIC - within_accX.ModelCriterion.AIC










% MOTION
selMot= isfinite(rtblMot.pt);


within_rt = fitlme(rtblMot, 'lrt ~ (targ+dist)  + (1+(targ+dist)  | pt)', 'Exclude', ~(selMot & rtblMot.acc == 1),...
    'FitMethod', 'REML', 'DummyVarCoding', 'effects');

[~,~,ffx_rt] = fixedEffects(within_rt, 'DFMethod', 'Satterthwaite')
[~,~,rfx_rt] = randomEffects(within_rt);
ffx_rt.D = round(ffx_rt.tStat ./ sqrt(ffx_rt.DF+1),3, 'significant');
ffx_rt.Lower = round(ffx_rt.Lower,3, 'significant');
ffx_rt.Upper = round(ffx_rt.Upper,3, 'significant')


within_acc = fitglme(rtblMot, 'acc ~ (targ+dist) + (1+(targ+dist) | pt)', 'Distribution', 'binomial', 'Exclude', ~(selMot & isfinite(rtblMot.rt)),...
    'FitMethod', 'Laplace', 'DummyVarCoding', 'effects');

[~,~,ffx_acc] = fixedEffects(within_acc);
ffx_acc.DF = repmat(npt-1, [3,1]);
ffx_acc.pValue = tcdf(abs(ffx_acc.tStat), npt-1, 'upper')*2
ffx_acc.D = round(ffx_acc.tStat ./ sqrt(npt),3, 'significant');
ffx_acc.Lower = round(ffx_acc.Lower,3, 'significant');
ffx_acc.Upper = round(ffx_acc.Upper,3, 'significant')

[~,~,rfx_acc] = randomEffects(within_acc);








%% adaptation

selCol = isfinite(rtbl.pt) & rtbl.acc1 == 1;


% no interaction
within_rt = fitlme(rtbl, 'lrt ~ targ + dist*dist1  + (1+ targ + dist*dist1  | pt)', 'Exclude', ~(selCol & rtbl.acc == 1),...
    'FitMethod', 'REML', 'DummyVarCoding', 'effects');

[~,~,ffx_rt] = fixedEffects(within_rt, 'DFMethod', 'Satterthwaite')


within_acc = fitglme(rtbl, 'acc ~ targ + dist*dist1 + (1 + targ + dist*dist1 | pt)', 'Distribution', 'binomial', 'Exclude', ~(selCol & isfinite(rtbl.rt)),...
    'FitMethod', 'Laplace', 'DummyVarCoding', 'effects');

[~,~,ffx_acc] = fixedEffects(within_acc);
ffx_acc.DF = repmat(npt-1, [3,1]);
ffx_acc.pValue = tcdf(abs(ffx_acc.tStat), npt-1, 'upper')*2




% with interaction
within_rtX = fitlme(rtbl, 'lrt ~ (targ+dist)*(targ1+dist1)  + (1+(targ*dist)*(targ1+dist1)  | pt)', 'Exclude', ~(selCol & rtbl.acc == 1),...
    'FitMethod', 'REML', 'DummyVarCoding', 'effects');
[~,~,ffx_rtX] = fixedEffects(within_rtX, 'DFMethod', 'Satterthwaite')


within_accX = fitglme(rtbl, 'acc ~ (targ*dist)*(targ1+dist1) + (1+(targ*dist)*(targ1+dist1) | pt)', 'Distribution', 'binomial', 'Exclude', ~(selCol & isfinite(rtbl.rt)),...
    'FitMethod', 'Laplace', 'DummyVarCoding', 'effects');


[~,~,ffx_accX] = fixedEffects(within_accX);
ffx_accX.DF = repmat(npt-1, [4,1]);
ffx_accX.pValue = tcdf(abs(ffx_accX.tStat), npt-1, 'upper')*2


% model comparision 
within_rt.ModelCriterion.AIC - within_rtX.ModelCriterion.AIC
within_acc.ModelCriterion.AIC - within_accX.ModelCriterion.AIC





%% within-trial dynamics

selCol = isfinite(rtbl.pt);


% no interaction
within_rt = fitlme(rtbl, 'lrt ~ (targ + dist)*acc  + (1+ (targ + dist)*acc  | pt)', 'Exclude', ~(selCol),...
    'FitMethod', 'REML', 'DummyVarCoding', 'effects');

[~,~,ffx_rt] = fixedEffects(within_rt, 'DFMethod', 'Satterthwaite')


within_acc = fitglme(rtbl, 'acc ~ (targ + dist)*crt + (1 + (targ + dist)*crt | pt)', 'Distribution', 'binomial', 'Exclude', ~(selCol & isfinite(rtbl.rt)),...
    'FitMethod', 'Laplace', 'DummyVarCoding', 'effects');

[~,~,ffx_acc] = fixedEffects(within_acc);
ffx_acc.DF = repmat(npt-1, [3,1]);
ffx_acc.pValue = tcdf(abs(ffx_acc.tStat), npt-1, 'upper')*2































%%



if length(unique(tcol.pt)) > 1

    % RT color
    mdl_color = fitlme(rtbl, 'lrt ~ (targ+dist)  + (1+(targ+dist)  | pt)', 'Exclude', ~(selCol & rtbl.acc == 1),...
        'FitMethod', 'REML', 'DummyVarCoding', 'effects')

    [~,~,rfx] = randomEffects(mdl_color);
    [ts,idx]=sort(rfx.tStat(ismember(rfx.Name, 'dist')))
    pts(idx)'

    % RT motion
    %     mdl_motion = fitlme(rtblMot, 'lrt ~ (targ+dist) + (1+(targ+dist) | pt)', 'Exclude', ~(selMot & rtblMot.acc == 1),...
    %         'FitMethod', 'REML', 'DummyVarCoding', 'effects')


    mdl = fitlme(rtbl, 'lrt ~ (targ+dist) * (targ1+dist1) + (1 + (targ+dist) * (targ1+dist1) | pt)', 'Exclude', ~(selCol & rtbl.acc == 1),...
        'FitMethod', 'REML', 'DummyVarCoding', 'effects')

    [~,~,rfx] = randomEffects(mdl);
    [ts,idx]=sort(rfx.tStat(ismember(rfx.Name, 'dist:dist1')))
    pts(idx)'


    % acc
    fitglme(rtbl, 'acc ~ (targ+dist) + (1+(targ+dist) | pt)', 'Distribution', 'binomial', 'Exclude', ~(selCol),...
        'FitMethod', 'Laplace', 'DummyVarCoding', 'effects')


    fitglme(rtbl, 'lapse ~ (targ+dist) + (1+(targ+dist) | pt)', 'Distribution', 'binomial', ...
        'FitMethod', 'Laplace', 'DummyVarCoding', 'effects')


    %     fitglme(rtblMot, 'acc ~ (targ+dist) + (1+(targ+dist) | pt)', 'Distribution', 'binomial', 'Exclude', ~(selMot),...
    %         'FitMethod', 'Laplace', 'DummyVarCoding', 'effects')

    fitglme(rtbl, 'acc ~ (targ+dist) * (targ1+dist1) + (1 + (targ+dist) * (targ1+dist1) | pt)', 'Distribution', 'binomial', 'Exclude', ~(selCol),...
        'FitMethod', 'Laplace', 'DummyVarCoding', 'effects')

else


    fitlm(rtbl, 'lrt ~ resp + (targ+dist)', 'Exclude', ~(selCol & rtbl.acc == 1), 'DummyVarCoding', 'effects')
    fitglm(rtbl, 'acc ~ resp + (targ+dist)', 'Distribution', 'binomial', 'Exclude', ~(selCol), 'DummyVarCoding', 'effects')


end




mdl_color = fitlme(rtbl, 'lrt ~ dur  + (1 +  dur  | pt)', 'Exclude', ~(selCol & rtbl.acc == 1),...
    'FitMethod', 'REML', 'DummyVarCoding', 'effects')

[~,~,rfx] = randomEffects(mdl_color);


figure;
plot(rfx.Estimate(ismember(rfx.Name, 'dur')) + rfx.Estimate(ismember(rfx.Name, '(Intercept)')), '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
yline(0, '-k', 'LineWidth', 1)
ylabel('RT ~ trial duration')
xlabel('participants')



%% PES


mdl_pes = fitlme(rtbl, 'lrt ~ pes + (1 + pes  | pt)', 'Exclude', ~(isfinite(rtbl.pt) & rtbl.acc == 1),...
    'FitMethod', 'REML', 'DummyVarCoding', 'effects')
[~,~,rfx] = randomEffects(mdl_pes);
figure; violinplot(mdl_pes.Coefficients.Estimate(ismember(mdl_pes.Coefficients.Name, 'pes')) + rfx.Estimate(ismember(rfx.Name, 'pes'))); yline(0)


fitglme(rtbl, 'acc ~ pes + (1 + rt*pes  | pt)', 'Distribution', 'binomial', 'Exclude', ~(isfinite(rtbl.pt)),...
    'FitMethod', 'Laplace', 'DummyVarCoding', 'effects')


%% plot combined
selCol = isfinite(tcol.pt);
selMot = isfinite(tmot.pt);



figure; hold on;
tiledlayout('flow', 'Padding', 'compact', 'TileSpacing', 'compact');


% ======== RT plot
nexttile; hold on;

[f,x] = ksdensity(tcol.rt(tcol.acc == 1 & selCol));
plot(x,f*mean(tcol.acc == 1), '-g', 'LineWidth', 2);
[f,x] = ksdensity(tcol.rt(tcol.acc == 1 & tcol.pt == pts(pp) & selCol));
plot(x,-(f*(1-mean(tcol.acc == 1))), '-r', 'LineWidth', 2);
yline(0, '-k', 'LineWidth', 2);
xline(1.5, '--k', 'LineWidth', 1);
xlabel('RT')
set(gca, 'TickDir', 'out', 'LineWidth', 1);


% ======== autocorr
nexttile; hold on;

yyaxis('left');
[f,x, bd] = autocorr(tcol.lrt(selCol), 5);
plot(x,f, '-or', 'markerfacecolor', 'r', 'LineWidth', 1, 'markersize', 6);
[f,x] = autocorr(tcol.acc,5);
plot(x,f, '-ob', 'markerfacecolor', 'b', 'LineWidth', 1, 'markersize', 6);

yline(0, '-k')
yline(bd(1), '--r')
yline(bd(2), '--r')


xlabel('lag')
ylabel('acf')
legend({'RT', 'Accuracy'}, 'location', 'northeast')
title('autocorr')


% ======== RT/accuracy over time
nexttile([1,2]); hold on;

yyaxis('right');
plot(movmean(tcol.rt,10,'omitnan'), '-r', 'LineWidth', 2)
ylabel('logRT (r)')

yyaxis('left');
plot(movmean(tcol.acc,10,'omitnan'), '-b', 'LineWidth', 2);
plot(movmean(isfinite(tcol.lrt),10,'omitnan'), '-m', 'LineWidth', 2);
ylabel('accuracy (b) / response (m)')

for ii = 1:6
    xline(150*ii)
end
xlim([0,length(tcol.acc)+1])
xlabel('trial')



% ======== Color Blocks

% targ
nexttile; hold on;

gt = grpstats(tcol(selCol,:), {'targ'}, {'mean', 'median'});

yyaxis('left');
plot(gt.targ, gt.mean_acc, '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 2);
ylabel('accuracy')
yyaxis('right');
plot(gt.targ, gt.median_lrt, '-or', 'MarkerFaceColor', 'r', 'LineWidth', 2);
ylabel('logRT')

xlim([.5, 5.5]);
xticks([])
xlabel('target coherence')
title('targ - attend-color')

set(gca, 'TickDir', 'out', 'LineWidth', 1);

% dist
nexttile; hold on;

gt = grpstats(tcol(selCol,:), {'dist'}, {'mean', 'median'});

yyaxis('left');
plot(gt.dist, gt.mean_acc, '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 2);
ylabel('accuracy')

yyaxis('right');
plot(gt.dist, gt.median_lrt, '-or', 'MarkerFaceColor', 'r', 'LineWidth', 2);
ylabel('logRT')

xlim([-2.5, 2.5]);
xticks([])
title('dist - attend-color')
xlabel('distractor congruence')

set(gca, 'TickDir', 'out', 'LineWidth', 1);


% ======== Motion Blocks

% targ
nexttile; hold on;

gt = grpstats(tmot(selMot,:), {'targ'}, {'mean', 'median'});

yyaxis('left');
plot(gt.targ, gt.mean_acc, '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 1);
yyaxis('right');
plot(gt.targ, gt.median_lrt, '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1);

title('targ - attend-motion')


% dist
nexttile; hold on;

gt = grpstats(tmot(selMot,:), {'dist'}, {'mean', 'median'});

yyaxis('left');
plot(gt.dist, gt.mean_acc, '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 1);
yyaxis('right');
plot(gt.dist, gt.median_lrt, '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1);

title('dist - attend-motion')




%% plot by run
selCol = isfinite(tcol.pt);

center = @(x) x-nanmean(x)

figure; hold on;

tcol.logRT = log(tcol.rt);

gt = grpstats(tcol(selCol,:), {'dist', 'run'}, {'mean', 'median'});

runs = unique(tcol.run);

cmap = colormap('parula');
cidx = round(linspace(1,256, length(runs)));


nexttile; hold on;
for ii = 1:length(runs)

    plot(gt.dist(gt.run == ii), center(gt.median_logRT(gt.run==ii)), '-', 'color', cmap(cidx(ii),:), 'LineWidth',2)

end
legend({'run1','run2','run3','run4','run5','run6'})


nexttile; hold on;
for ii = 1:length(runs)

    plot(gt.dist(gt.run == ii), center(gt.mean_acc(gt.run==ii)), '-', 'color', cmap(cidx(ii),:), 'LineWidth',2)

end
legend({'run1','run2','run3','run4','run5','run6'})







%% plot per subj

selCol = isfinite(tcol.pt);
selMot = isfinite(tmot.pt);


for pp = 1:npt


    figure; hold on;
    tiledlayout('flow', 'Padding', 'compact', 'TileSpacing', 'compact');


    % ======== RT plot
    nexttile; hold on;

    [f,x] = ksdensity(tcol.rt(tcol.acc == 1 & tcol.pt == pts(pp) & selCol));
    plot(x,f*mean(tcol.acc == 1), '-g', 'LineWidth', 2);
    [f,x] = ksdensity(tcol.rt(tcol.acc == 1 & tcol.pt == pts(pp) & selCol));
    plot(x,-(f*(1-mean(tcol.acc == 1))), '-r', 'LineWidth', 2);
    yline(0, '-k', 'LineWidth', 2);
    xline(1.5, '--k', 'LineWidth', 1);
    title(sprintf('pt %d', pp))
    xlabel('RT')


    % ======== autocorr
    nexttile; hold on;

    yyaxis('left');
    [f,x, bd] = autocorr(tcol.lrt(selCol), 5);
    plot(x,f, '-or', 'markerfacecolor', 'r', 'LineWidth', 1, 'markersize', 6);
    [f,x] = autocorr(tcol.acc,5);
    plot(x,f, '-ob', 'markerfacecolor', 'b', 'LineWidth', 1, 'markersize', 6);

    yline(0, '-k')
    yline(bd(1), '--r')
    yline(bd(2), '--r')


    xlabel('lag')
    ylabel('acf')
    legend({'RT', 'Accuracy'}, 'location', 'northeast')
    title('autocorr')


    % ======== RT/accuracy over time
    nexttile([1,2]); hold on;

    yyaxis('right');
    plot(movmean(tcol.rt(tcol.pt == pts(pp)),10,'omitnan'), '-r', 'LineWidth', 2)
    ylabel('logRT (r)')

    yyaxis('left');
    plot(movmean(tcol.acc(tcol.pt == pts(pp)),10,'omitnan'), '-b', 'LineWidth', 2);
    plot(movmean(isfinite(tcol.lrt(tcol.pt == pp)),10,'omitnan'), '-m', 'LineWidth', 2);
    ylabel('accuracy (b) / response (m)')

    for ii = 1:6
        xline(150*ii)
    end
    xlim([0,length(tcol.acc)+1])
    xlabel('trial')



    % ======== Color Blocks

    % targ
    nexttile; hold on;

    gt = grpstats(tcol(tcol.pt == pts(pp) & selCol,:), {'targ'}, {'mean', 'median'});

    yyaxis('left');
    plot(gt.targ, gt.mean_acc, '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 1);
    ylabel('accuracy')
    yyaxis('right');
    plot(gt.targ, gt.median_lrt, '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1);
    ylabel('logRT')

    title('targ - attend-color')


    % dist
    nexttile; hold on;

    gt = grpstats(tcol(tcol.pt == pts(pp) & selCol,:), {'dist'}, {'mean', 'median'});

    yyaxis('left');
    plot(gt.dist, gt.mean_acc, '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 1);
    yyaxis('right');
    plot(gt.dist, gt.median_lrt, '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1);

    title('dist - attend-color')



    % ======== Motion Blocks

    % targ
    nexttile; hold on;

    gt = grpstats(tmot(tmot.pt == pts(pp) & selMot,:), {'targ'}, {'mean', 'median'});

    yyaxis('left');
    plot(gt.targ, gt.mean_acc, '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 1);
    yyaxis('right');
    plot(gt.targ, gt.median_lrt, '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1);

    title('targ - attend-motion')


    % dist
    nexttile; hold on;

    gt = grpstats(tmot(tmot.pt == pts(pp) & selMot,:), {'dist'}, {'mean', 'median'});

    yyaxis('left');
    plot(gt.dist, gt.mean_acc, '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 1);
    yyaxis('right');
    plot(gt.dist, gt.median_lrt, '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1);

    title('dist - attend-motion')



    % stimulus confusion
    nexttile;
    stimMx = ones(4);
    gt = grpstats(tcol(tcol.pt == pts(pp) & selCol,:), {'stim', 'altStim'}, 'mean', 'DataVars', 'acc')
    for ii = 1:height(gt)
        stimMx(gt.stim(ii), gt.altStim(ii)) = gt.mean_acc(ii);
    end
    imagesc(stimMx'); colorbar;
    xlabel('colorTarget'); ylabel('colorAlt')
    xticks(1:4); xticklabels(d.data.results.rand.param.colorNames);
    yticks(1:4); yticklabels(d.data.results.rand.param.colorNames);
    title('color confusion')


    % stimulus sensitivity
    nexttile;
    gstim=grpstats(tcol(tcol.pt == pts(pp) & selCol,:), {'stim'}, 'mean', 'DataVars', 'acc');
    gAlt = grpstats(tcol(tcol.pt == pts(pp) & selCol,:), {'altStim'}, 'mean', 'DataVars', 'acc');
    bar(norminv(gstim.mean_acc) - (norminv(1-gAlt.mean_acc)));
    xticks(1:4); xticklabels(d.data.results.rand.param.colorNames);
    ylabel('d prime')
    title('color sensitivity')


end


























