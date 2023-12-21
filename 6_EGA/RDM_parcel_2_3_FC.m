%% analyze mediation
clear;clc;

%%
analysisFolder = 'parcelRSA_flipped' % 'parcelRSA_2022-06'; %'parcelRSA_orig'




analysisName = 'PFCl2IPS'




switch analysisName
    %% =================== KEY CONTRASTS

    case 'PFCl2IPS'
        %%


        directMdl1 = 'FC_PFCl_parcel'
        directMdl2 = 'FC_IPS_parcel'
        indirectMdl = 'FC_IPS-PFCl_parcel'



        conPairs1 = {...
            {'PFCl', 'cohTarg'},...
            {'PFCl', 'cohDist'},...
            {'PFCl', 'PFCl'},...
            {'cohTarg', 'cohTarg'},...
            {'cohDist', 'cohDist'},...
            };

        conNames1 = {...
            'PFCl_cohTarg',...
            'PFCl_cohDist',...
            'PFCl',...
            'cohTarg',...
            'cohDist',...
            };



        conPairs2 = {...
            {'IPS', 'cohTarg'},...
            {'IPS', 'cohDist'},...
            {'IPS', 'IPS'},...
            {'cohTarg', 'cohTarg'},...
            {'cohDist', 'cohDist'},...
            };

        conNames2 = {...
            'IPS_cohTarg',...
            'IPS_cohDist',...
            'IPS',...
            'cohTarg',...
            'cohDist',...
            };



        ddNames = {...
            'cohTarg',...
            'cohDist',...
            'PFCl',...
            'cohTargPure',...
            'cohDistPure',...
            };




end






% get paths
save_dir = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s',analysisFolder, analysisName)
data_dir1 = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s/fit-results',analysisFolder, directMdl1)
data_dir2 = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s/fit-results',analysisFolder, directMdl2)
data_ind = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s/fit-results',analysisFolder, indirectMdl)




mask = '/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii';
maskVo = spm_vol(mask);
maskImg = flipud(spm_read_vols(maskVo));
nParcel = 400;


% load data
dr_dir1 = dir(fullfile(data_dir1, '*.mat'));
dr_dir2 = dir(fullfile(data_dir2, '*.mat'));
dr_ind = dir(fullfile(data_ind, '*.mat'));

[R_dir1, R_dir2, R_ind] = deal([]);

for ii = 1:length(dr_dir1)

    f_dir1 = load(fullfile(dr_dir1(ii).folder, dr_dir1(ii).name));
    R_dir1 = cat(4, R_dir1, f_dir1.R);

    f_dir2 = load(fullfile(dr_dir2(ii).folder, dr_dir2(ii).name));
    R_dir2 = cat(4, R_dir2, f_dir2.R);


    f_ind = load(fullfile(dr_ind(ii).folder, dr_ind(ii).name));
    R_ind = cat(4, R_ind, f_ind.R);


end

conds_dir1 = f_dir1.Opt.condLabel
conds_dir2 = f_dir2.Opt.condLabel
conds_ind = f_ind.Opt.condLabel

nConds_dir1 = length(conds_dir1)
nConds_dir2 = length(conds_dir1)
nConds_ind = length(conds_ind)


npt = length(dr_dir1)


%% plot A-m-Y vs A-Y (main effect)
%
%
%
% alpha = .001;
% sel_dir = @(x) find(strcmp(conds_dir1, x));
% sel_ind = @(x) find(strcmp(conds_ind, x));
%
% mkdir(fullfile(save_dir, 'mediation_main'));
% delete(fullfile(save_dir, 'mediation_main', '*.nii'))
%
%
%
% % R reliability
% for vv = 1:length(mainNames)
%
%
%     [~,~,~,stats_dir] = ttest(squeeze(R_dir1(sel_dir(mainPairs{vv}{1}),sel_dir(mainPairs{vv}{1}),:,:))');
%     [~,~,~,stats_ind] = ttest(squeeze(R_ind(sel_ind(mainPairs{vv}{1}),sel_ind(mainPairs{vv}{1}),:,:))');
%
%     [~,~,~,stats_comp] = ttest(...
%         squeeze(R_ind(sel_ind(mainPairs{vv}{1}),sel_ind(mainPairs{vv}{1}),:,:))' - squeeze(R_dir1(sel_dir(mainPairs{vv}{1}),sel_dir(mainPairs{vv}{1}),:,:))'...
%         );
%
%
%     fname_dir = fullfile(save_dir, 'mediation_main', sprintf('dir_%s.nii', mainNames{vv}));
%     fname_ind = fullfile(save_dir, 'mediation_main', sprintf('ind_%s.nii', mainNames{vv}));
%     fname_comp = fullfile(save_dir, 'mediation_main', sprintf('comp_%s.nii', mainNames{vv}));
%
%
%     % ===== direct
%
%     % print to parcellation
%         parcelBrain = maskImg;
%     for rr = 1:nParcel
%
%         parcelBrain(parcelBrain == rr) = stats_dir.tstat(rr);
%
%     end
%
%
%
%     % write
%     VoOut      = struct(...
%         'fname',    fname_dir,...
%         'dim',      maskVo.dim,...
%         'dt',       [spm_type('float32') spm_platform('bigend')],...
%         'mat',      maskVo.mat,...
%         'n',        [1 1],...
%         'descrip',  'joint pattern reliability');
%
%     spm_write_vol(VoOut, parcelBrain);
%
%
%     % ===== ind
%
%     % print to parcellation
%         parcelBrain = maskImg;
%     for rr = 1:nParcel
%
%         parcelBrain(parcelBrain == rr) = stats_ind.tstat(rr);
%
%     end
%
%
%
%     % write
%     VoOut      = struct(...
%         'fname',    fname_ind,...
%         'dim',      maskVo.dim,...
%         'dt',       [spm_type('float32') spm_platform('bigend')],...
%         'mat',      maskVo.mat,...
%         'n',        [1 1],...
%         'descrip',  'joint pattern reliability');
%
%     spm_write_vol(VoOut, parcelBrain);
%
%
%
%
%     % ===== direct
%
%     % print to parcellation
%         parcelBrain = maskImg;
%     for rr = 1:nParcel
%
%         parcelBrain(parcelBrain == rr) = stats_comp.tstat(rr);
%
%     end
%
%
%
%     % write
%     VoOut      = struct(...
%         'fname',    fname_comp,...
%         'dim',      maskVo.dim,...
%         'dt',       [spm_type('float32') spm_platform('bigend')],...
%         'mat',      maskVo.mat,...
%         'n',        [1 1],...
%         'descrip',  'joint pattern reliability');
%
%     spm_write_vol(VoOut, parcelBrain);
%
%
%
% end



%% plot A-m-Y vs A-Y (alignment)

mx_sims = 10000;

hb_alpha = .05;
alpha = .001;

sel_dir1 = @(x) find(strcmp(conds_dir1, x));
sel_dir2 = @(x) find(strcmp(conds_dir2, x));
sel_ind = @(x) find(strcmp(conds_ind, x));

mkdir(fullfile(save_dir, 'mediation_con'));
delete(fullfile(save_dir, 'mediation_con', '*.nii'))



% R reliability
for vv = 1:length(conNames1)


    [~,~,~,stats_dir1] = ttest(squeeze(R_dir1(sel_dir1(conPairs1{vv}{1}),sel_dir1(conPairs1{vv}{2}),:,:))');
    [~,~,~,stats_dir2] = ttest(squeeze(R_dir2(sel_dir2(conPairs2{vv}{1}),sel_dir2(conPairs2{vv}{2}),:,:))');
    %     [~,~,~,stats_ind] = ttest(squeeze(R_ind(sel_ind(conPairs1{vv}{1}),sel_ind(conPairs1{vv}{2}),:,:))');

    dists_1 = ...
        squeeze(R_ind(sel_ind(conPairs1{vv}{1}),sel_ind(conPairs1{vv}{2}),:,:))' - ...
        squeeze(R_dir1(sel_dir1(conPairs1{vv}{1}),sel_dir1(conPairs1{vv}{2}),:,:))';
    [~,pval_comp1,~,stats_comp1] = ttest(dists_1);


    dists_2 = ...
        squeeze(R_ind(sel_ind(conPairs2{vv}{1}),sel_ind(conPairs2{vv}{2}),:,:))' - ...
        squeeze(R_dir2(sel_dir2(conPairs2{vv}{1}),sel_dir2(conPairs2{vv}{2}),:,:))';
    [~,pval_comp2,~,stats_comp2] = ttest(dists_2);


    dists_dd = ...
        (squeeze(R_ind(sel_ind(conPairs1{vv}{1}),sel_ind(conPairs1{vv}{2}),:,:))' - ...
        squeeze(R_dir1(sel_dir1(conPairs1{vv}{1}),sel_dir1(conPairs1{vv}{2}),:,:))') - ...
        ( squeeze(R_ind(sel_ind(conPairs2{vv}{1}),sel_ind(conPairs2{vv}{2}),:,:))' - ...
        squeeze(R_dir2(sel_dir2(conPairs2{vv}{1}),sel_dir2(conPairs2{vv}{2}),:,:))');


    [~,pval_dd,~,stats_dd] = ttest(dists_dd);




    fname_dir1 = fullfile(save_dir, 'mediation_con', sprintf('dir_%s.nii', conNames1{vv}));
    fname_dir2 = fullfile(save_dir, 'mediation_con', sprintf('dir_%s.nii', conNames2{vv}));

    %     fname_ind1 = fullfile(save_dir, 'mediation_con', sprintf('ind_%s.nii', conNames1{vv}));
    %     fname_ind2 = fullfile(save_dir, 'mediation_con', sprintf('ind_%s.nii', conNames2{vv}));

    fname_comp1 = fullfile(save_dir, 'mediation_con', sprintf('comp_%s.nii', conNames1{vv}));
    fname_comp1_hb = fullfile(save_dir, 'mediation_con', sprintf('comp_hb_%s.nii', conNames1{vv}));

    fname_comp2 = fullfile(save_dir, 'mediation_con', sprintf('comp_%s.nii', conNames2{vv}));
    fname_comp2_hb = fullfile(save_dir, 'mediation_con', sprintf('comp_hb_%s.nii', conNames2{vv}));

    fname_dd = fullfile(save_dir, 'mediation_con', sprintf('dd_%s.nii', ddNames{vv}));
    fname_dd_hb = fullfile(save_dir, 'mediation_con', sprintf('dd_hb_%s.nii', ddNames{vv}));





    % ===== direct 1 ======================================================

    % print to parcellation
    parcelBrain = maskImg;
    for rr = 1:nParcel

        parcelBrain(parcelBrain == rr) = stats_dir1.tstat(rr);

    end



    % write
    VoOut      = struct(...
        'fname',    fname_dir1,...
        'dim',      maskVo.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      maskVo.mat,...
        'n',        [1 1],...
        'descrip',  'joint pattern reliability');

    spm_write_vol(VoOut, parcelBrain);





    % ===== direct 2 ======================================================

    % print to parcellation
    parcelBrain = maskImg;
    for rr = 1:nParcel

        parcelBrain(parcelBrain == rr) = stats_dir2.tstat(rr);

    end



    % write
    VoOut      = struct(...
        'fname',    fname_dir2,...
        'dim',      maskVo.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      maskVo.mat,...
        'n',        [1 1],...
        'descrip',  'joint pattern reliability');

    spm_write_vol(VoOut, parcelBrain);







    % ===== contrast 1  ======================================================

    % print to parcellation
    parcelBrain = maskImg;
    for rr = 1:nParcel

        parcelBrain(parcelBrain == rr) = stats_comp1.tstat(rr);

    end



    % write
    VoOut      = struct(...
        'fname',    fname_comp1,...
        'dim',      maskVo.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      maskVo.mat,...
        'n',        [1 1],...
        'descrip',  'joint pattern reliability');

    spm_write_vol(VoOut, parcelBrain);





    % ===== contrast 1 (hb threshold)  ======================================================
    [~,pval1,~,~] = ttest(squeeze(R_ind(sel_ind(conPairs1{vv}{2}),sel_ind(conPairs1{vv}{2}),:,:))');
    pval2 = pval1;
    sel = pval1 <.001;

    % get correction
    %     pvalsel = (pval1<alpha) & (pval2<alpha);
    %     findPval = find(pvalsel);
    %     [psort, idxsort] = sort(pval_comp1(pvalsel));
    %     hb_pval = zeros(1, length(pval1));
    %     corP = psort .*[length(psort):-1:1];
    %     corP(find(corP>.05, 1):end) = 1;
    %     hb_pval(findPval(idxsort)) = corP;


    % randomization test
    perm_mx = reshape(datasample([-1,1], numel(dists_1)*mx_sims), [size(dists_1,1), size(dists_1,2), mx_sims]);
    null_R = repmat(dists_1, [1,1,mx_sims]) .* perm_mx;

    max_null_R = max(squeeze(mean(null_R(:,sel,:)) ./ std(null_R(:,sel,:))));
    min_null_R = min(squeeze(mean(null_R(:,sel,:)) ./ std(null_R(:,sel,:))));
    true_R = mean(dists_1)./std(dists_1);
    if isempty(max_null_R)
        max_null_R = 0;
        min_null_R = 0;
    end





    % print to parcellation
    parcelBrain = maskImg;
    for rr = 1:nParcel

        if true_R(rr) > 0
            parcelBrain(parcelBrain == rr) = (mean(true_R(rr)>max_null_R')>(1-(hb_alpha/2)));
        else
            parcelBrain(parcelBrain == rr) = (mean(true_R(rr)<min_null_R')>(1-(hb_alpha/2)));
        end
    end



    % write
    VoOut      = struct(...
        'fname',    fname_comp1_hb,...
        'dim',      maskVo.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      maskVo.mat,...
        'n',        [1 1],...
        'descrip',  'joint pattern reliability');

    spm_write_vol(VoOut, parcelBrain);










    % ===== contrast 2  ======================================================

    % print to parcellation
    parcelBrain = maskImg;
    for rr = 1:nParcel

        parcelBrain(parcelBrain == rr) = stats_comp2.tstat(rr);

    end



    % write
    VoOut      = struct(...
        'fname',    fname_comp2,...
        'dim',      maskVo.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      maskVo.mat,...
        'n',        [1 1],...
        'descrip',  'joint pattern reliability');

    spm_write_vol(VoOut, parcelBrain);





    % ===== contrast 2 (hb threshold)  ======================================================
    [~,pval1,~,~] = ttest(squeeze(R_ind(sel_ind(conPairs2{vv}{2}),sel_ind(conPairs2{vv}{2}),:,:))');
    pval2 = pval1;
    sel = pval1 <.001;

    % get correction
    %     pvalsel = (pval1<alpha) & (pval2<alpha);
    %     findPval = find(pvalsel);
    %     [psort, idxsort] = sort(pval_comp2(pvalsel));
    %     hb_pval = zeros(1, length(pval1));
    %     corP = psort .*[length(psort):-1:1];
    %     corP(find(corP>.05, 1):end) = 1;
    %     hb_pval(findPval(idxsort)) = corP;

    % randomization test
    perm_mx = reshape(datasample([-1,1], numel(dists_2)*mx_sims), [size(dists_2,1), size(dists_2,2), mx_sims]);
    null_R = repmat(dists_2, [1,1,mx_sims]) .* perm_mx;

    max_null_R = max(squeeze(mean(null_R(:,sel,:)) ./ std(null_R(:,sel,:))));
    min_null_R = min(squeeze(mean(null_R(:,sel,:)) ./ std(null_R(:,sel,:))));
    true_R = mean(dists_2)./std(dists_2);
    if isempty(max_null_R)
        max_null_R = 0;
        min_null_R = 0;
    end




    % print to parcellation
    parcelBrain = maskImg;
    for rr = 1:nParcel
        if true_R(rr) > 0
            parcelBrain(parcelBrain == rr) = (mean(true_R(rr)>max_null_R')>(1-(hb_alpha/2)));
        else
            parcelBrain(parcelBrain == rr) = (mean(true_R(rr)<min_null_R')>(1-(hb_alpha/2)));
        end
    end



    % write
    VoOut      = struct(...
        'fname',    fname_comp2_hb,...
        'dim',      maskVo.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      maskVo.mat,...
        'n',        [1 1],...
        'descrip',  'joint pattern reliability');

    spm_write_vol(VoOut, parcelBrain);








    % ===== DD  ======================================================

    % print to parcellation
    parcelBrain = maskImg;
    for rr = 1:nParcel

        parcelBrain(parcelBrain == rr) = stats_dd.tstat(rr);

    end



    % write
    VoOut      = struct(...
        'fname',    fname_dd,...
        'dim',      maskVo.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      maskVo.mat,...
        'n',        [1 1],...
        'descrip',  'joint pattern reliability');

    spm_write_vol(VoOut, parcelBrain);





    % ===== DD (hb threshold)  ======================================================
    [~,pval1,~,~] = ttest(squeeze(R_ind(sel_ind(conPairs1{vv}{2}),sel_ind(conPairs1{vv}{2}),:,:))');
    pval2 = pval1;
    sel = pval1 <.001;

    % get correction
    %     pvalsel = (pval1<alpha) & (pval2<alpha);
    %     findPval = find(pvalsel);
    %     [psort, idxsort] = sort(pval_dd(pvalsel));
    %     hb_pval = zeros(1, length(pval1));
    %     corP = psort .*[length(psort):-1:1];
    %     corP(find(corP>.05, 1):end) = 1;
    %     hb_pval(findPval(idxsort)) = corP;


    % randomization test
    perm_mx = reshape(datasample([-1,1], numel(dists_dd)*mx_sims), [size(dists_dd,1), size(dists_dd,2), mx_sims]);
    null_R = repmat(dists_dd, [1,1,mx_sims]) .* perm_mx;

    max_null_R = max(squeeze(mean(null_R(:,sel,:)) ./ std(null_R(:,sel,:))));
    min_null_R = min(squeeze(mean(null_R(:,sel,:)) ./ std(null_R(:,sel,:))));
    true_R = mean(dists_dd)./std(dists_dd);
    if isempty(max_null_R)
        max_null_R = 0;
        min_null_R = 0;
    end







    % print to parcellation
    parcelBrain = maskImg;
    for rr = 1:nParcel
        if true_R(rr) > 0
            parcelBrain(parcelBrain == rr) = (mean(true_R(rr)>max_null_R')>(1-(hb_alpha/2)));
        else
            parcelBrain(parcelBrain == rr) = (mean(true_R(rr)<min_null_R')>(1-(hb_alpha/2)));
        end
    end



    % write
    VoOut      = struct(...
        'fname',    fname_dd_hb,...
        'dim',      maskVo.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      maskVo.mat,...
        'n',        [1 1],...
        'descrip',  'joint pattern reliability');

    spm_write_vol(VoOut, parcelBrain);









end

fprintf('\n DONE CONTRAST \n')





%% bar plot mediation



sel_dir1 = @(x) find(strcmp(conds_dir1, x));
sel_dir2 = @(x) find(strcmp(conds_dir2, x));
sel_ind = @(x) find(strcmp(conds_ind, x));


roiList = {'mot', 'col', 'VisC', 'VisB', 'VisA'}
roiLoc = {[176, 375],...
    [181, 379],...
    [197, 198, 199, 200, 397, 398, 399, 400],...
    [185   186   187   188   189   190   191   192   193   194   195   196   384   385   386   387   388   389   390   391   392   393   394   395   396],...
    [   172   173   174   175   176   177   178   179   180   181   182   183   184   369   370   371   372   373   374   375   376   377   378   379   380   381   382   383],
    }


% roiList = {'mot', 'col', 'VisC'}
% roiLoc = {[176, 375],...
%     [181, 379],...
%     [197, 198, 199, 200, 397, 398, 399, 400],...
%     }



mkdir(fullfile(save_dir, 'plots'));
mkdir(fullfile(save_dir, 'regions'));
delete(fullfile(save_dir, 'regions', '*.nii'))




for rr = 1:length(roiList)

    regSel = ismember(1:400, roiLoc{rr});

    % make mask
    parcelBrain = zeros(size(maskImg));
    parcelBrain(ismember(maskImg, find(regSel))) = 1;


    fname = fullfile(save_dir, 'regions', sprintf('%s.nii', roiList{rr}));

    % write
    VoOut      = struct(...
        'fname',    fname,...
        'dim',      maskVo.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      maskVo.mat,...
        'n',        [1 1],...
        'descrip',  roiList{rr});

    spm_write_vol(VoOut, parcelBrain);

end




% R reliability
for vv = 1:length(conNames1)


    for rr = 1:length(roiList)


        regSel = ismember(1:400, roiLoc{rr});

        mRoi_dir1 = mean(mean(squeeze(R_dir1(sel_dir1(conPairs1{vv}{1}), sel_dir1(conPairs1{vv}{2}), regSel, :))));
        mRoi_dir2 = mean(mean(squeeze(R_dir2(sel_dir2(conPairs2{vv}{1}), sel_dir2(conPairs2{vv}{2}), regSel, :))));
        mRoi_ind1 = mean(mean(squeeze(R_ind(sel_ind(conPairs1{vv}{1}), sel_ind(conPairs1{vv}{2}), regSel, :))));
        mRoi_ind2 = mean(mean(squeeze(R_ind(sel_ind(conPairs2{vv}{1}), sel_ind(conPairs2{vv}{2}), regSel, :))));

        [~,p_1]  = ttest(...
            mean(squeeze(R_dir1(sel_dir1(conPairs1{vv}{1}),sel_dir1(conPairs1{vv}{2}),regSel,:))),...
            mean(squeeze(R_ind(sel_ind(conPairs1{vv}{1}),sel_ind(conPairs1{vv}{2}),regSel,:))));

        [~,p_2]  = ttest(...
            mean(squeeze(R_dir2(sel_dir2(conPairs2{vv}{1}),sel_dir2(conPairs2{vv}{2}),regSel,:))),...
            mean(squeeze(R_ind(sel_ind(conPairs2{vv}{1}),sel_ind(conPairs2{vv}{2}),regSel,:))));


        [~,p_x]  = ttest(...
            mean(squeeze(R_dir1(sel_dir1(conPairs1{vv}{1}),sel_dir1(conPairs1{vv}{2}),regSel,:))) - mean(squeeze(R_ind(sel_ind(conPairs1{vv}{1}),sel_ind(conPairs1{vv}{2}),regSel,:))),...
            mean(squeeze(R_dir2(sel_dir2(conPairs2{vv}{1}),sel_dir2(conPairs2{vv}{2}),regSel,:)))- mean(squeeze(R_ind(sel_ind(conPairs2{vv}{1}),sel_ind(conPairs2{vv}{2}),regSel,:))));



        sRoi_dir1 = std(mean(squeeze(R_dir1(sel_dir1(conPairs1{vv}{1}),sel_dir1(conPairs1{vv}{2}),regSel,:))))./sqrt(npt);
        sRoi_dir2 = std(mean(R_dir2(sel_dir2(conPairs2{vv}{1}),sel_dir2(conPairs2{vv}{2}),regSel,:)))./sqrt(npt);
        sRoi_ind1 = std(mean(R_ind(sel_ind(conPairs1{vv}{1}),sel_ind(conPairs1{vv}{2}),regSel,:)))./sqrt(npt);
        sRoi_ind2 = std(mean(R_ind(sel_ind(conPairs2{vv}{1}),sel_ind(conPairs2{vv}{2}),regSel,:)))./sqrt(npt);



        % plot
        f=figure; hold on;
        errorbar([.5, .6]-.01, [mRoi_dir1, mRoi_dir2], [sRoi_dir1, sRoi_dir2], 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2)
        errorbar([.5, .6]+.01, [mRoi_ind1, mRoi_ind2], [sRoi_ind1, sRoi_ind2], 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 15, 'LineWidth', 2)
        yline(0, '-k', 'LineWidth', 2)
        xlim([.45, .65])
        title([conNames1{vv}, '-', roiList{rr}])
        xticks([.5, .6]);
        set(gca, 'TickDir', 'out', 'LineWidth', 2)
        axis('square')
        ylabel('seed-coh corr')

        yl = ylim;
        try
            yticks([yl(1), 0, yl(2)])
        catch
            yticks([yl(1), yl(2)])
        end

        r1 = split(conNames1{1}, '_'); r2=split(conNames2{1}, '_');
        xticklabels({r1{1}, r2{1}})

        saveas(f, fullfile(save_dir, 'plots', [conNames1{vv}, '-', roiList{rr}]), 'pdf');


        fprintf('\n%s-%s_pval %f\n', [conNames1{vv}, '-', roiList{rr}],r1{1}, p_1)
        fprintf('%s-%s_pval %f\n\n', [conNames1{vv}, '-', roiList{rr}],r2{1}, p_2)
        fprintf('%s-X_pval %f\n\n', [conNames1{vv}, '-', roiList{rr}], p_x)


    end


end


close all;




%% test coherence and see correlation


