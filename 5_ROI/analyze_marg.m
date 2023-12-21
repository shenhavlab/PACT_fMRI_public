
clear; clc;

addpath(genpath('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util'))
addpath('/Users/hr0283/Documents/MATLAB/spm')

%% load data


load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/level 2/margCoh_s-8mm_cfd-movt_cvi-wls/kong22_dACC.mat');
[npt, nvx] = size(data{1})


 % mask
mask = load_nii('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/kong22_dACC.nii');
labels =  mask.img(mask.img>0);


ctrl_mask = load_nii('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/kongNets/kong_Cont.nii');
vattn_mask = load_nii('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/kongNets/kong_SalVenAttn.nii');


ctrl_idx =  ctrl_mask.img(mask.img>0);
vattn_idx =  vattn_mask.img(mask.img>0);


% make masks
vo = spm_vol('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/kongNets/kong_Cont.nii');
Vo      = struct(...
    'fname',    '/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/dACC_Cont.nii',...
    'dim',      vo.dim,...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      vo.mat,...
    'n',        [1 1],...
    'descrip',  'dACC Control');

spm_write_vol(Vo, ctrl_mask.img==1 & mask.img==1);


vo = spm_vol('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/kongNets/kong_SalVenAttn.nii');
Vo      = struct(...
    'fname',    '/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/dACC_SalVenAttn.nii',...
    'dim',      vo.dim,...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      vo.mat,...
    'n',        [1 1],...
    'descrip',  'dACC SalVenAttn');

spm_write_vol(Vo, vattn_mask.img==1 & mask.img==1);






% target sig
targSig = load_nii('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/level 2/feature_s-8mm_cfd-movt_cvi-wls/targCoh/dACC__tfce10000/T_tfce_tstat_fwep_c2.nii');
targSig = targSig.img(mask.img>0) >= 1.6;


% distractor sig
distSig = load_nii('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/level 2/feature_s-8mm_cfd-movt_cvi-wls/cong/dACC__tfce10000/T_tfce_tstat_fwep_c2.nii');
distSig = distSig.img(mask.img>0) >= 1.6;







%% exract and plot gradient

all_mds = cmdscale(pdist(mXYZ'), 2);


mds = all_mds(:,1);

xx = [];
for ii = 1:length(mds)
    xx(mXYZ(3,ii) - min(mXYZ(3,:)) + 1, mXYZ(2,ii) - min(mXYZ(2,:)) + 1) = mds(ii);

end

% figure;
% imagesc(flip(xx))
% axis('square');


%% make table
vec = @(x) x(:);
nanz = @(x,d) (x-nanmean(x,d))./nanstd(x,[],d);
cent = @(x) x-nanmean(x)

t = table;

t.pt = vec(repmat(1:npt, [nvx, 1]));
t.label = repmat(labels, [npt, 1]);
t.ctrl = repmat(ctrl_idx, [npt, 1]);
t.vsal = repmat(vattn_idx, [npt, 1]);




t.targSig = repmat(targSig, [npt, 1]);
t.distSig = repmat(distSig, [npt, 1]);


t.targ1 = vec(data{1}');
t.targ2 = vec(data{2}');
t.targ4 = vec(data{3}');
t.targ5 = vec(data{4}');

t.dist1 = vec(data{5}');
t.dist2 = vec(data{6}');
t.dist4 = vec(data{7}');
t.dist5 = vec(data{8}');




%% make plot


figure; hold on;
tiledlayout('flow', 'TileSpacing','compact', 'Padding','compact')

% full ROI

nexttile(1,[2,1]); hold on;
gt = grpstats(t, {'pt'});

allTarg = [gt.mean_targ1, gt.mean_targ2, gt.mean_targ4, gt.mean_targ5]; allTarg  = allTarg-mean(allTarg,2);
allDist = [gt.mean_dist1, gt.mean_dist2, gt.mean_dist4, gt.mean_dist5]; allDist  = allDist-mean(allDist,2);


errorbar([1,2,4,5]-.2, mean(allTarg), std(allTarg)./sqrt(npt), '-ok', 'MarkerFaceColor', 'auto', 'LineWidth',2)
errorbar([1,2,4,5]+.2, mean(allDist), std(allDist)./sqrt(npt), '-or', 'MarkerFaceColor', 'auto', 'LineWidth',2)
xlim([.5, 5.5])
title('full ROI')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
legend({'target', 'distractor'})


% % target vx
% nexttile; hold on;
% gt = grpstats(t(t.targSig==1,:), {'pt'});
% 
% allTarg = [gt.mean_targ1, gt.mean_targ2, gt.mean_targ4, gt.mean_targ5]; allTarg  = allTarg-mean(allTarg,2);
% allDist = [gt.mean_dist1, gt.mean_dist2, gt.mean_dist4, gt.mean_dist5]; allDist  = allDist-mean(allDist,2);
% 
% 
% plot([1,2,4,5]-.2, mean(allTarg),  '-ok', 'MarkerFaceColor', 'auto', 'LineWidth',2)
% errorbar([1,2,4,5]+.2, mean(allDist), std(allDist)./sqrt(npt), '-or', 'MarkerFaceColor', 'auto', 'LineWidth',2)
% xlim([.5, 5.5])
% title('target difficulty voxels')
% set(gca, 'TickDir', 'out', 'LineWidth', 1)
% 
% 
% 
% % distractor vx
% nexttile; hold on;
% gt = grpstats(t(t.distSig==1,:), {'pt'});
% 
% allTarg = [gt.mean_targ1, gt.mean_targ2, gt.mean_targ4, gt.mean_targ5]; allTarg  = allTarg-mean(allTarg,2);
% allDist = [gt.mean_dist1, gt.mean_dist2, gt.mean_dist4, gt.mean_dist5]; allDist  = allDist-mean(allDist,2);
% 
% 
% errorbar([1,2,4,5]-.2, mean(allTarg), std(allTarg)./sqrt(npt), '-ok', 'MarkerFaceColor', 'auto', 'LineWidth',2)
% plot([1,2,4,5]+.2, mean(allDist),  '-or', 'MarkerFaceColor', 'auto', 'LineWidth',2)
% xlim([.5, 5.5])
% title('distractor incongruence voxels')
% set(gca, 'TickDir', 'out', 'LineWidth', 1)



% outside control vx
nexttile; hold on;
gt = grpstats(t(t.vsal==1,:), {'pt'});

allTarg = [gt.mean_targ1, gt.mean_targ2, gt.mean_targ4, gt.mean_targ5]; allTarg  = allTarg-mean(allTarg,2);
allDist = [gt.mean_dist1, gt.mean_dist2, gt.mean_dist4, gt.mean_dist5]; allDist  = allDist-mean(allDist,2);


errorbar([1,2,4,5]-.2, mean(allTarg), std(allTarg)./sqrt(npt), '-ok', 'MarkerFaceColor', 'auto', 'LineWidth',2)
errorbar([1,2,4,5]+.2, mean(allDist), std(allDist)./sqrt(npt), '-or', 'MarkerFaceColor', 'auto', 'LineWidth',2)
xlim([.5, 5.5])
title('ventral attention network')
set(gca, 'TickDir', 'out', 'LineWidth', 1)


% within control vx
nexttile; hold on;
gt = grpstats(t(t.ctrl==1,:), {'pt'});

allTarg = [gt.mean_targ1, gt.mean_targ2, gt.mean_targ4, gt.mean_targ5]; allTarg  = allTarg-mean(allTarg,2);
allDist = [gt.mean_dist1, gt.mean_dist2, gt.mean_dist4, gt.mean_dist5]; allDist  = allDist-mean(allDist,2);


errorbar([1,2,4,5]-.2, mean(allTarg), std(allTarg)./sqrt(npt), '-ok', 'MarkerFaceColor', 'auto', 'LineWidth',2)
errorbar([1,2,4,5]+.2, mean(allDist), std(allDist)./sqrt(npt), '-or', 'MarkerFaceColor', 'auto', 'LineWidth',2)
xlim([.5, 5.5])
title('control network')
set(gca, 'TickDir', 'out', 'LineWidth', 1)




%% plot subj







figure; hold on;
tiledlayout('flow', 'TileSpacing','compact', 'Padding','compact')

% full ROI

nexttile(1,[2,1]); hold on;
gt = grpstats(t, {'pt'});

allTarg = [gt.mean_targ1, gt.mean_targ2, gt.mean_targ4, gt.mean_targ5]; allTarg  = allTarg-mean(allTarg,2);
allDist = [gt.mean_dist1, gt.mean_dist2, gt.mean_dist4, gt.mean_dist5]; allDist  = allDist-mean(allDist,2);


plot([1,2,4,5]-.2, (allTarg),  '-ok', 'MarkerFaceColor', 'auto', 'LineWidth',1.5)
plot([1,2,4,5]+.2, (allDist),  '-or', 'MarkerFaceColor', 'auto', 'LineWidth',1.5)
xlim([.5, 5.5])
title('full ROI')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
legend({'target', 'distractor'})


% target vx
nexttile; hold on;
gt = grpstats(t(t.targSig==1,:), {'pt'});

allTarg = [gt.mean_targ1, gt.mean_targ2, gt.mean_targ4, gt.mean_targ5]; allTarg  = allTarg-mean(allTarg,2);
allDist = [gt.mean_dist1, gt.mean_dist2, gt.mean_dist4, gt.mean_dist5]; allDist  = allDist-mean(allDist,2); 


plot([1,2,4,5]-.2, (allTarg),  '-ok', 'MarkerFaceColor', 'auto', 'LineWidth',1.5)
plot([1,2,4,5]+.2, (allDist),  '-or', 'MarkerFaceColor', 'auto', 'LineWidth',1.5)
xlim([.5, 5.5])
title('target difficulty voxels')
set(gca, 'TickDir', 'out', 'LineWidth', 1)



% distractor vx
nexttile; hold on;
gt = grpstats(t(t.distSig==1,:), {'pt'});

allTarg = [gt.mean_targ1, gt.mean_targ2, gt.mean_targ4, gt.mean_targ5]; allTarg  = allTarg-mean(allTarg,2);
allDist = [gt.mean_dist1, gt.mean_dist2, gt.mean_dist4, gt.mean_dist5]; allDist  = allDist-mean(allDist,2);


plot([1,2,4,5]-.2, (allTarg), '-ok', 'MarkerFaceColor', 'auto', 'LineWidth',1.5)
plot([1,2,4,5]+.2, (allDist),  '-or', 'MarkerFaceColor', 'auto', 'LineWidth',1.5)
xlim([.5, 5.5])
title('distractor incongruence voxels')
set(gca, 'TickDir', 'out', 'LineWidth', 1)













