%% RDM_parcel_2_1_analyze

clear;


%% SETUP analysis

analysisFolder = 'parcelRSA_flipped' % 'parcelRSA_2022-06'; %'parcelRSA_orig'

mask = '/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii';
maskVo = spm_vol(mask);
maskImg = flipud(spm_read_vols(maskVo));
nParcel = 400;




% get paths
analysisName = 'perfFullBlkTask_parcel' %


save_dir = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s',analysisFolder, analysisName)
data_dir = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s/fit-results',analysisFolder, analysisName)

dr = dir(fullfile(data_dir, '*.mat'));


% load data

[G_full,C_full,R_full,Rw] = deal([]);

for ii = 1:length(dr)

    f = load(fullfile(dr(ii).folder, dr(ii).name));

    %     if is_conn
    %
    %         R = cat(3, R, f.R);
    % %         Rw = cat(3, Rw, f.Rw);
    %
    %     else

    G_full = cat(4, G_full, f.G);
    C_full = cat(4, C_full, f.C);
    R_full = cat(4, R_full, f.R);

    %     end



end

condLabel = f.Opt.condLabel
nCond = length(condLabel)
npt = length(dr)



% get paths
% analysisFolder = 'parcelRSA_2022-06' % 'parcelRSA_2022-06'; %'parcelRSA_orig'
% analysisName = 'perfFullBlkTask_parcel' %


analysisName = 'perfBlkTask_parcel' %


save_dir = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s',analysisFolder, analysisName)
data_dir = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s/fit-results',analysisFolder, analysisName)

dr = dir(fullfile(data_dir, '*.mat'));


% load data

[G_match,C_match,R_match,Rw] = deal([]);

for ii = 1:length(dr)

    f = load(fullfile(dr(ii).folder, dr(ii).name));

    %     if is_conn
    %
    %         R = cat(3, R, f.R);
    % %         Rw = cat(3, Rw, f.Rw);
    %
    %     else

    G_match = cat(4, G_match, f.G);
    C_match = cat(4, C_match, f.C);
    R_match = cat(4, R_match, f.R);

    %     end



end

condLabel = f.Opt.condLabel
nCond = length(condLabel)
npt = length(dr)









%% =========




accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];
roiList = {'dACC', 'PFCl', 'SPL', 'IPS'}

condSel = @(x) find(strcmp(condLabel, x));

taskCols = {[20, 70, 160]/255, [219, 48, 105]/255};

for nn = 1

    figure;
    tiledlayout('flow', 'padding', 'compact', 'TileSpacing','compact')

    switch nn

        case 1

            regionList = roiList;
            regionTbl = f.info.region;

        case 2

            regionList = subnetworkList;
            regionTbl = f.info.subNetwork;

        case 3

            regionList = networkList;
            regionTbl = f.info.network;

    end


    for rr = 1:length(regionList)

        switch regionList{rr}
            case 'dACC'
                regSel = ismember(1:400, accSel);
            otherwise
                regSel = cellfun(@any, regexpi(regionTbl, regionList{rr}));
        end


        col_ms = [...
            squeeze(mean(R_full(condSel('colRespTarg'),condSel('colRespTarg'),regSel,:),3)),...
            squeeze(mean(R_full(condSel('colRespDist'),condSel('colRespDist'),regSel,:),3)),...
            squeeze(mean(R_full(condSel('colCohTarg'),condSel('colCohTarg'),regSel,:),3)),...
            squeeze(mean(R_full(condSel('colCohDist'),condSel('colCohDist'),regSel,:),3)),...
            squeeze(mean(R_full(condSel('colCong'),condSel('colCong'),regSel,:),3)),...
            ];


         col_match_ms = [...
            squeeze(mean(R_match(condSel('colRespTarg'),condSel('colRespTarg'),regSel,:),3)),...
            squeeze(mean(R_match(condSel('colRespDist'),condSel('colRespDist'),regSel,:),3)),...
            squeeze(mean(R_match(condSel('colCohTarg'),condSel('colCohTarg'),regSel,:),3)),...
            squeeze(mean(R_match(condSel('colCohDist'),condSel('colCohDist'),regSel,:),3)),...
            squeeze(mean(R_match(condSel('colCong'),condSel('colCong'),regSel,:),3)),...
            ];



        mot_ms = [...
            squeeze(mean(R_full(condSel('motRespTarg'),condSel('motRespTarg'),regSel,:),3)),...
            squeeze(mean(R_full(condSel('motRespDist'),condSel('motRespDist'),regSel,:),3)),...
            squeeze(mean(R_full(condSel('motCohTarg'),condSel('motCohTarg'),regSel,:),3)),...
            squeeze(mean(R_full(condSel('motCohDist'),condSel('motCohDist'),regSel,:),3)),...
            squeeze(mean(R_full(condSel('motCong'),condSel('motCong'),regSel,:),3)),...
            ];


        disp(regionList{rr})
        [~,col_pval,~,col_stats] = ttest(col_ms)
        [~,mot_pval,~,mot_stats] = ttest(mot_ms)
        [~,comp_pval,~,comp_stats] = ttest(col_ms, mot_ms)
        [~,resp_pval,~,resp_stats] = ttest((col_ms(:,1) + mot_ms(:,1)), (col_ms(:,2) + mot_ms(:,2)))


        nexttile; hold on;
        errorbar([1:5]-.2, mean(col_ms), std(col_ms)./sqrt(npt),...
            'o', 'Color', taskCols{1}, 'MarkerFaceColor', taskCols{1}, 'MarkerSize', 8, 'LineWidth', 1.5);

        errorbar([1:5]-.1, mean(col_match_ms), std(col_match_ms)./sqrt(npt),...
            'o', 'Color', taskCols{1}, 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5);

        errorbar([1:5]+.2, mean(mot_ms), std(mot_ms)./sqrt(npt),...
            'o', 'Color', taskCols{2}, 'MarkerFaceColor', taskCols{2},'MarkerSize', 8, 'LineWidth', 1.5);

        title(regionList{rr})
        xlim([.5, 5.5])
        xticks(1:5)
        xticklabels(condLabel(1:end-1))
        yline(0)
        ylim([-.02, .1])
        set(gca, 'TickDir', 'out', 'LineWidth', 1)





    end






end





