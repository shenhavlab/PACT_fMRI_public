%% analyze gradients

clear; clc;

addpath(genpath('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util'))


%% load data

load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/level 2/feature_s-8mm_cfd-movt_cvi-wls/kong22_dACC.mat');


ctrlNii = load_nii('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/kongNets/kong_Cont.nii');
ctrlCNii = load_nii('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/kongNets/kong_ContC.nii');


mask = load_nii('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/kong22_dACC.nii');
maskImg = mask.img;


ctrl = double(ctrlNii.img(maskImg>0));
ctrlC = double(ctrlCNii.img(maskImg>0));



[npt, nvx] = size(data{1})





%% extract and plot gradient

% all_mds = cmdscale(pdist(mXYZ'), 2);

[l,s,r] = pca(mXYZ(2:3,:)');
all_mds = s;

mds = all_mds(:,1);

xx = [];
for ii = 1:length(mds)
    xx(mXYZ(3,ii) - min(mXYZ(3,:)) + 1, mXYZ(2,ii) - min(mXYZ(2,:)) + 1) = mds(ii);

end

figure;
imagesc(flip(xx))
axis('square');


%% make table
vec = @(x) x(:);
center = @(x) x - nanmean(x);

nDisc = 10;

t = table;
t.mX =  repmat(mXYZ(1,:)', [npt, 1]);
t.mY =  repmat(mXYZ(2,:)', [npt, 1]);
t.mZ =  repmat(mXYZ(3,:)', [npt, 1]);
t.mds = repmat(mds, [npt, 1]);
t.mds2 = repmat(all_mds(:,2), [npt, 1]);
t.dmds = repmat(discretize(mds,nDisc), [npt, 1]);
t.targ = vec(data{3}');
t.dist = vec(data{5}');
t.distCoh = vec(data{4}');
t.cTarg = vec(center(data{3}'));
t.cDist = vec(center(data{5}'));
t.pt = vec(repmat(1:npt, [nvx, 1]));
% t.label = repmat(labels, [npt, 1]);
t.ctrl = repmat(center(ctrl), [npt, 1]);
t.ctrlC = repmat(center(ctrlC), [npt, 1]);
t.signX = t.mX > 45; t.signX(t.signX==0) = -1;

% t = t(t.mX > 45, :);


clear rDiff rCoh
for pp = 1:size(data{3},1)

    rDiff(pp) = corr(data{3}(pp,:)', data{5}(pp,:)', 'rows', 'pairwise');
    rCoh(pp) = corr(data{3}(pp,:)', data{4}(pp,:)', 'rows', 'pairwise');

end


log10(bf_ttest(rDiff))
log10(bf_ttest(rCoh))




%% do regression

targM   = fitlme(t, 'targ ~ mds^2 + (1 + mds^2 |pt)', 'FitMethod', 'REML');
[~,~,ffx_targ] = fixedEffects(targM, 'DFMethod', 'satterthwaite')
[~,~,rfx_targ] = randomEffects(targM);

ffx_targ.D = round(ffx_targ.tStat ./ sqrt(ffx_targ.DF+1),3, 'significant');
ffx_targ.Lower = round(ffx_targ.Lower,3, 'significant')
ffx_targ.Upper = round(ffx_targ.Upper,3, 'significant')

distM   = fitlme(t, 'dist ~ mds^2 + (1 + mds^2  |pt)', 'FitMethod', 'REML');
[~,~,ffx_dist] = fixedEffects(distM, 'DFMethod', 'satterthwaite')
[~,~,rfx_dist] = randomEffects(distM);


ffx_dist.D = round(ffx_dist.tStat ./ sqrt(ffx_dist.DF+1),3, 'significant');
ffx_dist.Lower = round(ffx_dist.Lower,3, 'significant')
ffx_dist.Upper = round(ffx_dist.Upper,3, 'significant')


tgM     = fitlme(t, 'targ ~ dist  + (1 + dist  |pt)', 'FitMethod', 'REML');
[~,~,ffx_cross] = fixedEffects(tgM, 'DFMethod', 'satterthwaite')
ffx_cross.D = round(ffx_cross.tStat ./ sqrt(ffx_cross.DF+1),3, 'significant');
ffx_cross.Lower = round(ffx_cross.Lower,3, 'significant')
ffx_cross.Upper = round(ffx_cross.Upper,3, 'significant')

log10(bf_ttest('T',ffx_cross.tStat(end), 'N', 29))

t.predTarg = targM.predict;
t.predDist = distM.predict;


%% regression discont
opt = statset('LinearMixedModel');
opt.UseParallel = true;


% parametric 
tic
targM0   = fitlme(t, 'targ ~ mds^2 + (1 + mds^2 |pt)', 'FitMethod', 'ML');
distM0   = fitlme(t, 'dist ~ mds^2 + (1 + mds^2  |pt)', 'FitMethod', 'ML');
toc;

targM1  = fitlme(t, 'targ ~ ctrlC + mds^2 + (1 + ctrlC + mds^2 |pt)', 'FitMethod', 'ML');
[~,~,ffx_targC] = fixedEffects(targM1, 'DFMethod', 'satterthwaite')
[table_targNest, siminfo_targNest] = compare(targM0, targM1, 'CheckNesting', true)
diff(table_targNest.BIC)
% tic
% [table_targ, siminfo_targ] = compare(targM0, targM1,'nsim',1000,'Options',opt, 'CheckNesting', true)
% toc

%     Simulated Likelihood Ratio Test: Nsim = 1000, Alpha = 0.05
% 
%     Model     DF    AIC            BIC            LogLik        LRStat    pValue      Lower         Upper    
%     targM0    10    -4.0642e+05    -4.0632e+05    2.0322e+05                                                 
%     targM1    15    -4.0773e+05    -4.0758e+05    2.0388e+05    1320.5    0.000999    2.5292e-05    0.0055534




distM1   = fitlme(t, 'dist ~ ctrlC + mds^2 + (1 + ctrlC + mds^2  |pt)', 'FitMethod', 'ML');
[~,~,ffx_distC] = fixedEffects(distM1, 'DFMethod', 'satterthwaite')
[table_distNest, siminfo_distNest] = compare(distM0, distM1, 'CheckNesting', true)
diff(table_distNest.BIC)

% tic
% [table_dist, siminfo_dist] = compare(distM0, distM1,'nsim',1000,'Options',opt, 'CheckNesting', true)
% toc


%     Simulated Likelihood Ratio Test: Nsim = 1000, Alpha = 0.05
% 
%     Model     DF    AIC            BIC            LogLik        LRStat    pValue      Lower         Upper    
%     distM0    10    -5.7113e+05    -5.7103e+05    2.8557e+05                                                 
%     distM1    15    -5.7427e+05    -5.7412e+05    2.8715e+05    3148.8    0.000999    2.5292e-05    0.0055534




% vs control-only 
tic
targM0   = fitlme(t, 'targ ~ ctrlC + (1 + ctrlC |pt)', 'FitMethod', 'ML');
distM0   = fitlme(t, 'dist ~ ctrlC + (1 + ctrlC  |pt)', 'FitMethod', 'ML');
toc;

targM1  = fitlme(t, 'targ ~ ctrlC + mds^2 + (1 + ctrlC + mds^2 |pt)', 'FitMethod', 'ML');
[~,~,ffx_targC] = fixedEffects(targM1, 'DFMethod', 'satterthwaite')
[table_targNest, siminfo_targNest] = compare(targM0, targM1, 'CheckNesting', true)
diff(table_targNest.BIC)

% tic
% [table_targ, siminfo_targ] = compare(targM0, targM1,'nsim',1000,'Options',opt, 'CheckNesting', true)
% toc

%     Simulated Likelihood Ratio Test: Nsim = 1000, Alpha = 0.05
% 
%     Model     DF    AIC            BIC            LogLik        LRStat    pValue      Lower         Upper    
%     targM0    10    -4.0642e+05    -4.0632e+05    2.0322e+05                                                 
%     targM1    15    -4.0773e+05    -4.0758e+05    2.0388e+05    1320.5    0.000999    2.5292e-05    0.0055534




distM1   = fitlme(t, 'dist ~ ctrlC + mds^2 + (1 + ctrlC + mds^2  |pt)', 'FitMethod', 'ML');
[~,~,ffx_distC] = fixedEffects(distM1, 'DFMethod', 'satterthwaite')
[table_distNest, siminfo_distNest] = compare(distM0, distM1, 'CheckNesting', true)
diff(table_distNest.BIC)

% tic
% [table_dist, siminfo_dist] = compare(distM0, distM1,'nsim',1000,'Options',opt, 'CheckNesting', true)
% toc


%     Simulated Likelihood Ratio Test: Nsim = 1000, Alpha = 0.05
% 
%     Model     DF    AIC            BIC            LogLik        LRStat    pValue      Lower         Upper    
%     distM0    10    -5.7113e+05    -5.7103e+05    2.8557e+05                                                 
%     distM1    15    -5.7427e+05    -5.7412e+05    2.8715e+05    3148.8    0.000999    2.5292e-05    0.0055534








%% full model
compM   = fitlme(t, 'mds ~ (cTarg^2 + cDist^2) + (1 + (cTarg^2 + cDist^2) | pt)', 'FitMethod', 'REML');
[~,~,ffx_comp] = fixedEffects(compM, 'DFMethod', 'satterthwaite')
[a,b,rfx_comp] = randomEffects(compM);


ffx_comp.D = round(ffx_comp.tStat ./ sqrt(ffx_comp.DF+1),3, 'significant');
ffx_comp.Lower = round(ffx_comp.Lower,3, 'significant')
ffx_comp.Upper = round(ffx_comp.Upper,3, 'significant')


t.predMDS = compM.predict;
for ii = 1:npt
t.predMDS(t.pt == ii) = discretize(t.predMDS(t.pt == ii),nDisc);
end


%% get minima

getB_ffx = @(name) ffx_comp.Estimate(ismember(ffx_comp.Name, name));
getB = @(name) ffx_comp.Estimate(ismember(ffx_comp.Name, name)) + rfx_comp.Estimate(ismember(rfx_comp.Name, name))

%  y = ax^2 + bx + c, min = c - b^2/4a.

targ_min = getB('(Intercept)') - (getB('cTarg').^2)./(4*getB('cTarg^2'));
dist_min = getB('(Intercept)') - (getB('cDist').^2)./(4*getB('cDist^2'));

targ_min_ffx = getB_ffx('(Intercept)') - (getB_ffx('cTarg').^2)./(4*getB_ffx('cTarg^2'));
dist_min_ffx = getB_ffx('(Intercept)') - (getB_ffx('cDist').^2)./(4*getB_ffx('cDist^2'));

targMinRank = 1+mean(t.mds < mean(targ_min))*(nDisc-1)
distMinRank = 1+mean(t.mds < mean(dist_min))*(nDisc-1)


% [~,idx] = min(abs(mean(targ_min_ffx) - t.mds));
% targMinRank = t.dmds(idx);
% 
% [~,idx] = min(abs(mean(dist_min_ffx) - t.mds));
% distMinRank = t.dmds(idx);


[~,p,~,stats] = ttest(targ_min - dist_min)
[p,~,rankstats]=signrank(targ_min - dist_min)
[rval,pval] = corr(targ_min, dist_min)
log10(bf_corr(rval, 29))

figure; hold on;
ksdensity(targ_min-dist_min)
xline(0)


%% make plot




% plot gradient

figure;
nexttile;
imagesc(flip(xx))
% axis('square');




% plot effect

nexttile; hold on;


gt = grpstats(t, {'dmds', 'pt'});
gt = grpstats(gt, {'dmds'}, {'mean', 'sem'});

% gtp = grpstats(t, {'predMDS', 'pt'});
% gtp = grpstats(gtp, {'predMDS'}, {'mean', 'sem'});


errorbar(gt.dmds + .1, gt.mean_mean_targ, gt.sem_mean_targ,...
    'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1)

errorbar(gt.dmds - .1, gt.mean_mean_dist, gt.sem_mean_dist,...
    'or', 'MarkerFaceColor', 'r', 'LineWidth', 1)


plot(gt.dmds + .1, gt.mean_mean_predTarg,...
    '-k', 'LineWidth', 2)

plot(gt.dmds - .1, gt.mean_mean_predDist,...
    '-r', 'LineWidth', 2)

% plot(gtp.predMDS + .1, gtp.mean_mean_targ,...
%     '-k', 'LineWidth', 2)
% 
% plot(gtp.predMDS - .1, gtp.mean_mean_dist,...
%     '-r', 'LineWidth', 2)


xline(targMinRank + .1, '-k', 'LineWidth', 1)
xline(distMinRank - .1, '-r', 'LineWidth', 1)


yline(0);

legend({'target', 'distractor'}, 'Location', 'northwest')

xlim([min(gt.dmds)-.5, max(gt.dmds)+.5])

set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlabel('gradient quantile')
ylabel('coeffecient')






    %% analyze breakpoint





%     % plot img
%     axis('equal')
%     title(contrast)
%     hold on;
% 
% 
%     % plot Sal boundary
%     [~,f]=contour(ismember(maskImg, [311, 312, 314, 107, 108, 110]), 1, '-w');
%     f.LineWidth = 1;


