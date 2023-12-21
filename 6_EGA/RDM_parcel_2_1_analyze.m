%% RDM_parcel_2_1_analyze

clear;
conPairs = [];

addpath(genpath('/Users/hr0283/Documents/MATLAB/spm12'))
addpath(genpath('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util'))


% SETUP analysis
analysisName = 'perf_parcel' %

analysisFolder = 'parcelRSA_flipped' % 'parcelRSA_2022-06'; %'parcelRSA_orig'



switch analysisName
 case 'FC_PFCl_parcel'
        %%
        % correlations
        corrPairs = {...
            {'respTarg', 'PFCl'},...
            {'cohTarg', 'PFCl'},...
            {'cohDist', 'PFCl'},...
            {'cong', 'PFCl'},...
            }

        corrName = {'respTarg', 'targCoh', 'distCoh', 'cong'}



    case {'FC_IPS-PFCl_parcel'}
        %%

        % correlations
        corrPairs = {...
            {'cohTarg', 'IPS'},...
            {'cohDist', 'IPS'},...
            {'cohTarg', 'PFCl'},...
            {'cohDist', 'PFCl'},...
            {'IPS', 'PFCl'},...
            }

        corrName = {...
            'targCohIPS', 'distCohIPS', ...
            'targCohdPFCl', 'distCohdPFCl', ...
            'IPS-PFCl', ...
            }



    case {'FC_IPS_parcel'}
        %%

        % correlations
        corrPairs = {...
            {'respTarg', 'IPS'},...
            {'respDist', 'IPS'},...
            {'cohTarg', 'IPS'},...
            {'cohDist', 'IPS'},...
            {'cong', 'IPS'},...
            }

        corrName = {'respTarg','respDist','targCoh', 'distCoh', 'cong'}



   case {'perfFullBlkTask_parcel'}
        %%

        % correlations
        corrPairs = {...
            {'colRespTarg', 'motRespTarg'},...
            {'colRespDist', 'motRespDist'},...
            {'colCohTarg', 'motCohTarg'},...
            {'colCohDist', 'motCohDist'},...
            {'colCong', 'motCong'},...
            {'colCohTarg', 'motCohDist'},...
            {'colCohDist', 'motCohTarg'},...
            }

        corrName = {'respTarg', 'respDist', 'cohTarg', 'cohDist', 'cong', 'color', 'motion'}

    case {'margCohTask_parcel'}
        %%
        corrPairs = {...
            {'targ4L', 'targ4R'},...
            {'dist2L', 'dist2R'},...
            }

        corrName = {'targResp', 'distResp'}

    case {'feature_parcel_conn'}
        %%

        % correlations
        corrPairs = {...
            {'cohTarg', 'conn1'},...
            {'cohTarg', 'conn2'},...
            {'cohDist', 'conn1'},...
            {'cohDist', 'conn2'},...
            {'cong', 'conn1'},...
            {'cong', 'conn2'},...
            }

        corrName = {...
            'targ1', 'targ2',...
            'dist1', 'dist2',...
            'cong1', 'cong2',...
            }


   case {'perfBlkTask_parcel'}
        %%

        % correlations
        corrPairs = {...
            {'colRespTarg', 'motRespTarg'},...
            {'colRespDist', 'motRespDist'},...
            {'colCohTarg', 'motCohTarg'},...
            {'colCohDist', 'motCohDist'},...
            {'colCong', 'motCong'},...
            {'colCohTarg', 'motCohDist'},...
            {'colCohDist', 'motCohTarg'},...
            {'colRespTarg', 'motRespDist'},...
            {'colRespDist', 'motRespTarg'},...
            }

        corrName = {'respTarg', 'respDist', 'cohTarg', 'cohDist', 'cong', 'colorCoh', 'motionCoh', 'colorResp', 'motionResp'}




   case {'featureBlkTask_parcel'}
        %%

        % correlations
        corrPairs = {...
            {'colRespTarg', 'motRespTarg'},...
            {'colRespDist', 'motRespDist'},...
            {'colCohTarg', 'motCohTarg'},...
            {'colCohDist', 'motCohDist'},...
            {'colCong', 'motCong'},...
            {'colCohTarg', 'motCohDist'},...
            {'colCohDist', 'motCohTarg'},...
            }

        corrName = {'respTarg', 'respDist', 'cohTarg', 'cohDist', 'cong', 'color', 'motion'}





    case {'margResp_parcel'}
        %%
        corrPairs = {...
            {'targ4L', 'targ4R'},...
            {'dist2L', 'dist2R'},...
            }

        corrName = {'targResp', 'distResp'}


    case {'feature_parcel'}
        %%

        % correlations
        corrPairs = {...
            {'respTarg', 'respDist'},...
            {'cohTarg', 'cohDist'},...
            {'respTarg', 'cohTarg'},...
            {'respDist', 'cohDist'},...
            {'respDist', 'cong'},...
            {'cohDist', 'cong'},...
            {'cohTarg', 'cong'},...
            }

        corrName = {...
            'resp', 'coh',...
            'targ', 'dist',...
            'respCong', 'cohCong',...
            'targCong',...
            }




    case {'perf_parcel'}
        %%

        % correlations
        corrPairs = {...
            {'respTarg', 'rt'},...
            {'respDist', 'rt'},...
            {'respTarg', 'acc'},...
            {'respDist', 'acc'},...
            {'cohTarg', 'rt'},...
            {'cohDist', 'rt'},...
            {'cohTarg', 'acc'},...
            {'cohDist', 'acc'},...
            {'cong', 'rt'},...
            {'cong', 'acc'},...
            {'rt', 'acc'}...
            {'cohTarg', 'cohDist'}...
            {'respTarg', 'respDist'}...
            {'respTarg', 'cohTarg'}...
            }

        corrName = {...
            'respTargRt', 'respDistRt',...
            'respTargAcc', 'respDistAcc',...
            'cohTargRt', 'cohDistRt',...
            'cohTargAcc', 'cohDistAcc',...
            'congRt', 'congAcc',...
            'rtAcc',...
            'coh', 'resp',...
            'repCohTarg'
            }



  end


% get paths
save_dir = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s',analysisFolder, analysisName)
data_dir = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s/fit-results',analysisFolder, analysisName)

dr = dir(fullfile(data_dir, '*.mat'));



mask = '/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/masks/Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii';
maskVo = spm_vol(mask);
maskImg = flipud(spm_read_vols(maskVo)); % fuck, need to flip parcel
nParcel = 400;


% load data
is_conn = any(regexp(analysisName, 'conn'));

[G,C,R,T,S,absR] = deal([]);

for ii = 1:length(dr)

    f = load(fullfile(dr(ii).folder, dr(ii).name));

    G = cat(4, G, f.G);
    C = cat(4, C, f.C);
    R = cat(4, R, f.R);
    try
        T = cat(4, T, f.T);
        S = cat(4, S, f.S);
        absR = cat(4, absR, f.absR);
    catch
    end

    %     end



end

condLabel = f.Opt.condLabel
nCond = length(condLabel)
npt = length(dr)


%% make overall reliability images

hb_alpha = .05;
mx_sims = 1e4;
mkdir(fullfile(save_dir, 'patternRel'));
delete(fullfile(save_dir, 'patternRel', '*.nii'))


for mm = 3 % for each metric

    % G reliability
    for vv = 1:size(G,1)

        % get stats
        switch mm

            case 1
                dists = squeeze(G(vv,vv,:,:))';
                mR = nanmean(squeeze(G(vv,vv,:,:))');
                [~,~,~,stats] = ttest(squeeze(G(vv,vv,:,:))');
                fname = fullfile(save_dir,'patternRel', sprintf('relG_%s.nii', condLabel{vv}));
                fnameR = fullfile(save_dir,'patternRel', sprintf('relG-stat_%s.nii', condLabel{vv}));
                fnameMX = fullfile(save_dir,'patternRel', sprintf('mx__relG_%s.nii', condLabel{vv}));
                fnameBF = fullfile(save_dir,'patternRel', sprintf('relG-bf_%s.nii', condLabel{vv}));


            case 2
                dists = squeeze(C(vv,vv,:,:))';
                mR = nanmean(squeeze(C(vv,vv,:,:))');
                [~,~,~,stats] = ttest(squeeze(C(vv,vv,:,:))');
                fname = fullfile(save_dir,'patternRel', sprintf('relC_%s.nii', condLabel{vv}));
                fnameR = fullfile(save_dir,'patternRel', sprintf('relC-stat_%s.nii', condLabel{vv}));
                fnameMX = fullfile(save_dir,'patternRel', sprintf('mx__relC_%s.nii', condLabel{vv}));
                fnameBF = fullfile(save_dir,'patternRel', sprintf('relC-bf_%s.nii', condLabel{vv}));


            case 3
                dists = squeeze(R(vv,vv,:,:))';
                mR = nanmean(squeeze(R(vv,vv,:,:))');
                [~,~,~,stats] = ttest(squeeze(R(vv,vv,:,:))');
                fname = fullfile(save_dir,'patternRel', sprintf('relR_%s.nii', condLabel{vv}));
                fnameR = fullfile(save_dir,'patternRel', sprintf('relR-stat_%s.nii', condLabel{vv}));
                fnameMX = fullfile(save_dir,'patternRel', sprintf('mx__relR_%s.nii', condLabel{vv}));
                fnameBF = fullfile(save_dir,'patternRel', sprintf('relR-bf_%s.nii', condLabel{vv}));


        end


        % reliablity T
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = stats.tstat(rr);
        end

        % write
        VoOut      = struct(...
            'fname',    fname,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  'pattern reliability');

        spm_write_vol(VoOut, parcelBrain);





        % reliablity R stat
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = mR(rr);
        end

        % write
        VoOut      = struct(...
            'fname',    fnameR,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  'pattern reliability');

        spm_write_vol(VoOut, parcelBrain);




        % reliablity BF
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = log10(bf_ttest(dists(:, rr)));
        end

        % write
        VoOut      = struct(...
            'fname',    fnameBF,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  'pattern reliability');

        spm_write_vol(VoOut, parcelBrain);







        % randomization test
        perm_mx = reshape(datasample([-1,1], numel(dists)*mx_sims), [size(dists,1), size(dists,2), mx_sims]);
        null_R = repmat(dists, [1,1,mx_sims]) .* perm_mx;

        max_null_R = max(squeeze(mean(null_R) ./ std(null_R)));
        min_null_R = min(squeeze(mean(null_R) ./ std(null_R)));
        true_R = mean(dists)./std(dists);
        if isempty(max_null_R)
            max_null_R = 0;
            min_null_R = 0;
        end


        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            if true_R(rr) > 0
                parcelBrain(parcelBrain == rr) = (mean(true_R(rr)>max_null_R')>(1-(hb_alpha/2)));
            else
                parcelBrain(parcelBrain == rr) = (mean(true_R(rr)<min_null_R')>(1-(hb_alpha/2)));
            end
        end

        % write
        VoOut      = struct(...
            'fname',    fnameMX,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  'pattern reliability');

        spm_write_vol(VoOut, parcelBrain);








    end
end

fprintf('\ndone: reliability\n')




%% make joint reliability images

alpha = .001;
condSel = @(x) find(strcmp(condLabel, x));

mkdir(fullfile(save_dir, 'jointRel'));
delete(fullfile(save_dir, 'jointRel', '*.nii'))


for mm = [3] % for each metric

    % G reliability
    for vv = 1:length(corrName)

        % get stats
        switch mm
            case 1

                [~,pval1,~,~] = ttest(squeeze(G(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{1}),:,:))');
                [~,pval2,~,~] = ttest(squeeze(G(condSel(corrPairs{vv}{2}),condSel(corrPairs{vv}{2}),:,:))');

                fname = fullfile(save_dir,'jointRel', sprintf('jointRelG_%s.nii', corrName{vv}));

            case 2

                [~,pval1,~,~] = ttest(squeeze(C(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{1}),:,:))');
                [~,pval2,~,~] = ttest(squeeze(C(condSel(corrPairs{vv}{2}),condSel(corrPairs{vv}{2}),:,:))');
                fname = fullfile(save_dir,'jointRel', sprintf('jointRelC_%s.nii', corrName{vv}));

            case 3

                [~,pval1,~,~] = ttest(squeeze(R(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{1}),:,:))');
                [~,pval2,~,~] = ttest(squeeze(R(condSel(corrPairs{vv}{2}),condSel(corrPairs{vv}{2}),:,:))');
                fname = fullfile(save_dir,'jointRel', sprintf('jointRelR_%s.nii', corrName{vv}));

        end


        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = (pval1(rr)<alpha) & (pval2(rr)<alpha);
        end

        % write
        VoOut      = struct(...
            'fname',    fname,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  'joint pattern reliability');

        spm_write_vol(VoOut, parcelBrain);



        % get join parcels
        if strcmp(corrName{vv}, 'coh') & mm == 3

            [pmin ] = min(pval1(cellfun(@any, regexpi(f.info.region, 'SPL'))) .* pval2(cellfun(@any, regexpi(f.info.region, 'SPL'))));
            SPL_coh_max = find(pval1.*pval2 == pmin & cellfun(@any, regexpi(f.info.region, 'SPL'))')
            [pmin] = min(pval1(cellfun(@any, regexpi(f.info.region, 'IPS'))) .* pval2(cellfun(@any, regexpi(f.info.region, 'IPS'))));
            IPS_coh_max = find(pval1.*pval2 == pmin & cellfun(@any, regexpi(f.info.region, 'IPS'))')


            SPL_coh = find(pval1'<alpha & pval2'<alpha & cellfun(@any, regexpi(f.info.region, 'SPL')))'
            IPS_coh = find(pval1'<alpha & pval2'<alpha & cellfun(@any, regexpi(f.info.region, 'IPS')))'
            PFCl_coh = find(pval1'<alpha & pval2'<alpha & cellfun(@any, regexpi(f.info.region, 'PFCl')))'

        end

    end
end


fprintf('\ndone: joint reliability\n')




%% correlate patterns

condSel = @(x) find(strcmp(condLabel, x));

mkdir(fullfile(save_dir, 'corrs'));
delete(fullfile(save_dir, 'corrs', '*.nii'))

clear all_R all_dR
for mm = 3 % for each metric

    % G reliability
    for vv = 1:length(corrName)

        % get stats
        switch mm

            %             case 1
            %
            %                 [~,~,~,stats] = ttest(squeeze(G(condSel(jointRelPairs{vv}{1}),condSel(jointRelPairs{vv}{2}),:,:))');
            %                 m = mean(squeeze(G(condSel(jointRelPairs{vv}{1}),condSel(jointRelPairs{vv}{2}),:,:))');
            %
            %                 fname = fullfile(save_dir,'corrs', sprintf('G_%s.nii', jointRelName{vv}));
            %
            %             case 2
            %
            %                 [~,pval1,~,~] = ttest(squeeze(C(condSel(jointRelPairs{vv}{1}),condSel(jointRelPairs{vv}{1}),:,:))');
            %                 [~,pval2,~,~] = ttest(squeeze(C(condSel(jointRelPairs{vv}{2}),condSel(jointRelPairs{vv}{2}),:,:))');
            %                 fname = fullfile(save_dir,'jointRel', sprintf('jointRelC_%s.nii', jointRelName{vv}));

            case 3

                % original corr
                con_R = squeeze(R(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{2}),:,:))';
                [~,~,ci,stats] = ttest(con_R);
                fname_M = fullfile(save_dir, 'corrs', sprintf('rM_%s.nii', corrName{vv}));
                fname_T = fullfile(save_dir, 'corrs', sprintf('rT_%s.nii', corrName{vv}));
                fname_BF = fullfile(save_dir, 'corrs', sprintf('rBF_%s.nii', corrName{vv}));

                all_R(:,vv) = mean(con_R);


                % deattenuated corr
                rel1 = squeeze(R(condSel(corrPairs{vv}{1}), condSel(corrPairs{vv}{1}),:,:));
                rel2 = squeeze(R(condSel(corrPairs{vv}{2}), condSel(corrPairs{vv}{2}),:,:));
                rel12 = squeeze(R(condSel(corrPairs{vv}{1}), condSel(corrPairs{vv}{2}),:,:));

                con_dR = (rel12./sqrt(rel1.*rel2))';
                con_dR(imag(con_dR)~= 0) = nan;

                [~,~,ci_d,stats_d] = ttest(con_dR);


                fname_dM = fullfile(save_dir, 'corrs', sprintf('drM_%s.nii', corrName{vv}));
                fname_dT = fullfile(save_dir, 'corrs', sprintf('drT_%s.nii', corrName{vv}));
                fname_dBF = fullfile(save_dir, 'corrs', sprintf('drBF_%s.nii', corrName{vv}));


                all_dR(:,vv) = mean(con_dR);


        end



        % =====  M (orig) ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = nanmean(con_R(:,rr));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_M,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);


        % =====  M (de-atten)  ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = nanmean(con_dR(:,rr));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_dM,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);




        % =====  T stat (orig) ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = stats.tstat(rr);
        end

        % write
        VoOut      = struct(...
            'fname',    fname_T,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);


        % =====  T stat (de-atten) ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = stats_d.tstat(rr);
        end

        % write
        VoOut      = struct(...
            'fname',    fname_dT,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);




        % =====  BF (orig) ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = log10(bf_ttest(con_R(:,rr)));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_BF,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);


        % =====  BF (de-atten)  ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = log10(bf_ttest(con_dR(:,rr)));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_dBF,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);




    end
end


fprintf('\ndone: correlated patterns\n')



%% correlate patterns (masked by reliability)


alpha = .001;
condSel = @(x) find(strcmp(condLabel, x));

mkdir(fullfile(save_dir, 'maskCorrs'));
delete(fullfile(save_dir, 'maskCorrs', '*.nii'))

clear all_R all_dR
for mm = 1:3 % for each metric

    % G reliability
    for vv = 1:length(corrName)

        % get stats
        switch mm

            case 1

                dists = S;

                fname_M = fullfile(save_dir, 'maskCorrs', sprintf('spearM_%s.nii', corrName{vv}));
                fname_T = fullfile(save_dir, 'maskCorrs', sprintf('spearT_%s.nii', corrName{vv}));
                fname_BF = fullfile(save_dir, 'maskCorrs', sprintf('spearBF_%s.nii', corrName{vv}));

                fname_dM = fullfile(save_dir, 'maskCorrs', sprintf('dspearM_%s.nii', corrName{vv}));
                fname_dT = fullfile(save_dir, 'maskCorrs', sprintf('dspearT_%s.nii', corrName{vv}));
                fname_dBF = fullfile(save_dir, 'maskCorrs', sprintf('dspearBF_%s.nii', corrName{vv}));
                fname_rel = fullfile(save_dir, 'maskCorrs', sprintf('dspearRel_%s.nii', corrName{vv}));
                fname_N = fullfile(save_dir, 'maskCorrs', sprintf('dspearN_%s.nii', corrName{vv}));


            case 2

                dists = T;

                fname_M = fullfile(save_dir, 'maskCorrs', sprintf('cosM_%s.nii', corrName{vv}));
                fname_T = fullfile(save_dir, 'maskCorrs', sprintf('cosT_%s.nii', corrName{vv}));
                fname_BF = fullfile(save_dir, 'maskCorrs', sprintf('cosBF_%s.nii', corrName{vv}));

                fname_dM = fullfile(save_dir, 'maskCorrs', sprintf('dcosM_%s.nii', corrName{vv}));
                fname_dT = fullfile(save_dir, 'maskCorrs', sprintf('dcosT_%s.nii', corrName{vv}));
                fname_dBF = fullfile(save_dir, 'maskCorrs', sprintf('dcosBF_%s.nii', corrName{vv}));
                fname_rel = fullfile(save_dir, 'maskCorrs', sprintf('dcosRel_%s.nii', corrName{vv}));
                fname_N = fullfile(save_dir, 'maskCorrs', sprintf('dcosN_%s.nii', corrName{vv}));

            case 3

                dists = R;

                fname_M = fullfile(save_dir, 'maskCorrs', sprintf('rM_%s.nii', corrName{vv}));
                fname_T = fullfile(save_dir, 'maskCorrs', sprintf('rT_%s.nii', corrName{vv}));
                fname_BF = fullfile(save_dir, 'maskCorrs', sprintf('rBF_%s.nii', corrName{vv}));

                fname_dM = fullfile(save_dir, 'maskCorrs', sprintf('drM_%s.nii', corrName{vv}));
                fname_dT = fullfile(save_dir, 'maskCorrs', sprintf('drT_%s.nii', corrName{vv}));
                fname_dBF = fullfile(save_dir, 'maskCorrs', sprintf('drBF_%s.nii', corrName{vv}));
                fname_rel = fullfile(save_dir, 'maskCorrs', sprintf('drRel_%s.nii', corrName{vv}));
                fname_N = fullfile(save_dir, 'maskCorrs', sprintf('drN_%s.nii', corrName{vv}));

        end


        % joint reliability mask
        [~,pval1,~,~] = ttest(squeeze(dists(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{1}),:,:))');
        [~,pval2,~,~] = ttest(squeeze(dists(condSel(corrPairs{vv}{2}),condSel(corrPairs{vv}{2}),:,:))');


        % original corr
        con_R = squeeze(dists(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{2}),:,:))';
        [~,~,ci,stats] = ttest(con_R);


        all_R(:,vv) = mean(con_R);


        % deattenuated corr
        rel1 = squeeze(dists(condSel(corrPairs{vv}{1}), condSel(corrPairs{vv}{1}),:,:));
        rel2 = squeeze(dists(condSel(corrPairs{vv}{2}), condSel(corrPairs{vv}{2}),:,:));
        rel12 = squeeze(dists(condSel(corrPairs{vv}{1}), condSel(corrPairs{vv}{2}),:,:));

        %                 jrel = max(1e-3^2, rel1.*rel2);
        jrel = rel1.*rel2;
        jrel(jrel<=0) = nan;
        con_dR = (rel12./sqrt(jrel))';
        con_dR(imag(con_dR)~= 0) = nan;

        [~,~,ci_d,stats_d] = ttest(con_dR);




        all_dR(:,vv) = mean(con_dR);



        % =====  M (orig) ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = nanmean(con_R(:,rr)) .* ((pval1(rr)<alpha) & (pval2(rr)<alpha));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_M,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);


        % =====  M (de-atten)  ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = nanmean(con_dR(:,rr)) .* ((pval1(rr)<alpha) & (pval2(rr)<alpha));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_dM,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);




        % =====  T stat (orig) ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = stats.tstat(rr) .* ((pval1(rr)<alpha) & (pval2(rr)<alpha));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_T,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);


        % =====  T stat (de-atten) ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = stats_d.tstat(rr) .* ((pval1(rr)<alpha) & (pval2(rr)<alpha));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_dT,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);




        % =====  BF (orig) ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = log10(bf_ttest(con_R(:,rr))) .* ((pval1(rr)<alpha) & (pval2(rr)<alpha));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_BF,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);


        % =====  BF (de-atten)  ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = log10(bf_ttest(con_dR(:,rr))) .* ((pval1(rr)<alpha) & (pval2(rr)<alpha));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_dBF,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);




        % =====  reliability  ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = nanmedian(sqrt(jrel(rr,:)')) .* ((pval1(rr)<alpha) & (pval2(rr)<alpha));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_rel,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);



        % =====  count  ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = mean(isfinite(con_dR(:,rr))) .* ((pval1(rr)<alpha) & (pval2(rr)<alpha));
        end

        % write
        VoOut      = struct(...
            'fname',    fname_N,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);




    end
end



fprintf('\ndone: correlated patterns (Masked by reliability)\n')







%% correlate patterns (masked by reliability & HB corrected)


mx_sims = 10000;

hb_alpha = .05;
alpha = .001;
condSel = @(x) find(strcmp(condLabel, x));

mkdir(fullfile(save_dir, 'threshMaskCorrs'));
delete(fullfile(save_dir, 'threshMaskCorrs', '*.nii'))

clear all_R all_dR
for mm = 2:3 % for each metric

    % G reliability
    for vv = 1:length(corrName)

        % get stats
        switch mm

            case 1

                dists = G;

                fname_T = fullfile(save_dir, 'threshMaskCorrs', sprintf('hb__gT_%s.nii', corrName{vv}));
                fname_bin = fullfile(save_dir, 'threshMaskCorrs', sprintf('hb__gBin_%s.nii', corrName{vv}));

            case 2

                dists = C;

                fname_T = fullfile(save_dir, 'threshMaskCorrs', sprintf('hb__cT_%s.nii', corrName{vv}));
                fname_bin = fullfile(save_dir, 'threshMaskCorrs', sprintf('hb__cBin_%s.nii', corrName{vv}));

                fname_mxT = fullfile(save_dir, 'threshMaskCorrs', sprintf('mx__cT_%s.nii', corrName{vv}));
                fname_mxbin = fullfile(save_dir, 'threshMaskCorrs', sprintf('mx__cBin_%s.nii', corrName{vv}));


            case 3

                dists = R;



                fname_T = fullfile(save_dir, 'threshMaskCorrs', sprintf('hb__rT_%s.nii', corrName{vv}));
                fname_bin = fullfile(save_dir, 'threshMaskCorrs', sprintf('hb__rBin_%s.nii', corrName{vv}));


                fname_mxT = fullfile(save_dir, 'threshMaskCorrs', sprintf('mx__rT_%s.nii', corrName{vv}));
                fname_mxbin = fullfile(save_dir, 'threshMaskCorrs', sprintf('mx__rBin_%s.nii', corrName{vv}));



        end



        % joint reliability mask
        [~,pval1,~,~] = ttest(squeeze(dists(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{1}),:,:))');
        [~,pval2,~,~] = ttest(squeeze(dists(condSel(corrPairs{vv}{2}),condSel(corrPairs{vv}{2}),:,:))');


        % original corr
        con_R = squeeze(dists(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{2}),:,:))';
        [~,cor_pval,~,stats] = ttest(con_R);

        % get correction
        pvalsel = pval1<alpha & pval2<alpha;
        findPval = find(pvalsel);
        [psort, idxsort] = sort(cor_pval(pvalsel));
        hb_pval = zeros(1, length(pval1));
        corP = psort .*[length(psort):-1:1];
        corP(find(corP>hb_alpha, 1):end) = 1;
        hb_pval(findPval(idxsort)) = corP;


        % randomization test
        perm_mx = reshape(datasample([-1,1], numel(con_R(:,pvalsel))*mx_sims), [size(con_R,1), sum(pvalsel), mx_sims]);
        null_R = repmat(con_R(:,pvalsel), [1,1,mx_sims]) .* perm_mx;

        max_null_R = max(squeeze(mean(null_R) ./ std(null_R)));
        min_null_R = min(squeeze(mean(null_R) ./ std(null_R)));
        true_R = mean(con_R)./std(con_R);
        if isempty(max_null_R)
            max_null_R = 0;
            min_null_R = 0;
        end









        % =====  T stat (orig) ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = stats.tstat(rr) .* ((pval1(rr)<alpha) & (pval2(rr)<alpha)) .* (hb_pval(rr)< hb_alpha);
        end

        % write
        VoOut      = struct(...
            'fname',    fname_T,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);




        % =====  T bin (orig) ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            parcelBrain(parcelBrain == rr) = ((pval1(rr)<alpha) & (pval2(rr)<alpha)) .* (hb_pval(rr)< hb_alpha);
        end

        % write
        VoOut      = struct(...
            'fname',    fname_bin,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);






        % =====  T bin (orig) ============================================
        parcelBrain = maskImg;

        % print to parcellation
        for rr = 1:nParcel
            if true_R(rr) > 0
                parcelBrain(parcelBrain == rr) = ((pval1(rr)<alpha) & (pval2(rr)<alpha)) .* (mean(true_R(rr)>max_null_R')>(1-(hb_alpha/2)));
            else
                parcelBrain(parcelBrain == rr) = ((pval1(rr)<alpha) & (pval2(rr)<alpha)) .* (mean(true_R(rr)<min_null_R')>(1-(hb_alpha/2)));
            end
        end

        % write
        VoOut      = struct(...
            'fname',    fname_mxbin,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  corrName{vv});

        spm_write_vol(VoOut, parcelBrain);






    end
end





fprintf('\n\nFINISHED MAIN ANALYSES\n\n')







%% correlate patterns (different levels of prior)


if strcmp(analysisName, 'feature_parcel')


    scales = logspace(log10(.707*.5), log10(.707*2), 5)

    alpha = .001;
    condSel = @(x) find(strcmp(condLabel, x));

    mkdir(fullfile(save_dir, 'bfScale'));
    delete(fullfile(save_dir, 'bfScale', '*.nii'))


    clear all_R all_dR
    for mm = 3 % for each metric

        % G reliability
        for vv = 1:length(corrName)

            % get stats
            switch mm

                %             case 1
                %
                %                 [~,~,~,stats] = ttest(squeeze(G(condSel(jointRelPairs{vv}{1}),condSel(jointRelPairs{vv}{2}),:,:))');
                %                 m = mean(squeeze(G(condSel(jointRelPairs{vv}{1}),condSel(jointRelPairs{vv}{2}),:,:))');
                %
                %                 fname = fullfile(save_dir,'corrs', sprintf('G_%s.nii', jointRelName{vv}));
                %
                %             case 2
                %
                %                 [~,pval1,~,~] = ttest(squeeze(C(condSel(jointRelPairs{vv}{1}),condSel(jointRelPairs{vv}{1}),:,:))');
                %                 [~,pval2,~,~] = ttest(squeeze(C(condSel(jointRelPairs{vv}{2}),condSel(jointRelPairs{vv}{2}),:,:))');
                %                 fname = fullfile(save_dir,'jointRel', sprintf('jointRelC_%s.nii', jointRelName{vv}));

                case 3

                    dists = R;

                    % joint reliability mask
                    [~,pval1,~,~] = ttest(squeeze(dists(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{1}),:,:))');
                    [~,pval2,~,~] = ttest(squeeze(dists(condSel(corrPairs{vv}{2}),condSel(corrPairs{vv}{2}),:,:))');


                    % original corr
                    con_R = squeeze(dists(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{2}),:,:))';

            end


            for pr = 1:length(scales)

                fname_BF = fullfile(save_dir, 'bfScale', sprintf('rBF_%s_%s.nii', corrName{vv}, num2str(scales(pr),2)));


                % =====  BF (orig) ============================================
                parcelBrain = maskImg;

                % print to parcellation
                for rr = 1:nParcel
                    parcelBrain(parcelBrain == rr) = log10(bf_ttest(con_R(:,rr), 'scale', scales(pr))) .* ((pval1(rr)<alpha) & (pval2(rr)<alpha));
                end

                % write
                VoOut      = struct(...
                    'fname',    fname_BF,...
                    'dim',      maskVo.dim,...
                    'dt',       [spm_type('float32') spm_platform('bigend')],...
                    'mat',      maskVo.mat .* [-1 1 1 1]',...
                    'n',        [1 1],...
                    'descrip',  corrName{vv});

                spm_write_vol(VoOut, parcelBrain);

            end




        end
    end



end



fprintf('\ndone: correlated patterns (Masked by reliability + FWE)\n')









%% contrast images


if ~isempty(conPairs)

    alpha = .001;
    condSel = @(x) find(strcmp(condLabel, x));

    mkdir(fullfile(save_dir, 'contrasts'));
    delete(fullfile(save_dir, 'contrasts', '*.nii'))


    for mm = [3] % for each metric

        % G reliability
        for vv = 1:length(conName)

            % get stats
            switch mm
                case 1

                    [~,pval1,~,~] = ttest(squeeze(G(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{1}),:,:))');
                    [~,pval2,~,~] = ttest(squeeze(G(condSel(corrPairs{vv}{2}),condSel(corrPairs{vv}{2}),:,:))');

                    fname = fullfile(save_dir,'contrasts', sprintf('conG_%s.nii', corrName{vv}));

                case 2

                    [~,pval1,~,~] = ttest(squeeze(C(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{1}),:,:))');
                    [~,pval2,~,~] = ttest(squeeze(C(condSel(corrPairs{vv}{2}),condSel(corrPairs{vv}{2}),:,:))');
                    fname = fullfile(save_dir,'contrasts', sprintf('conC_%s.nii', corrName{vv}));

                case 3

                    [~,pval1,~,stats] = ttest(squeeze(R(condSel(conPairs{vv}{1}),condSel(conPairs{vv}{1}),:,:))' - squeeze(R(condSel(conPairs{vv}{2}),condSel(conPairs{vv}{2}),:,:))');
                    fname = fullfile(save_dir,'contrasts', sprintf('conR_%s.nii', conName{vv}));

            end


            parcelBrain = maskImg;

            % print to parcellation
            for rr = 1:nParcel
                parcelBrain(parcelBrain == rr) = stats.tstat(rr);
            end

            % write
            VoOut      = struct(...
                'fname',    fname,...
                'dim',      maskVo.dim,...
                'dt',       [spm_type('float32') spm_platform('bigend')],...
                'mat',      maskVo.mat .* [-1 1 1 1]',...
                'n',        [1 1],...
                'descrip',  'contrast');

            spm_write_vol(VoOut, parcelBrain);



        end
    end



end









%% plot by region
accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];
fefSel = [315, 251, 116, 58];

PFClSalSel = [103, 104, 105, 307, 308, 309];
PFClCtrlSel = [59, 60, 61, 62, 252, 253, 254, 268, 278, 279, 280];
PFClCtrlSel = [278   279   280];
% PFClCtrlSel = [59, 60, 61, 62, 252, 253, 254, 268];


accSalSel = [    85    86   106   286   287   295];
accCtrlSel = [    73    80    81   255   281   282   283];


motSel = [176, 375];
visaSel =  [172   173   174   175   176   177   178   179   180   181   182   369   370   371   372   373   374   375   376   377   378   379   380];
M1LeftSel = [149, 150, 151];
PrCdSel = [   116   130   131   315   329   330   331];
PFCdSel = [     5    20    21    22    23    32    33    45    46    58    72    92   115   203   217   218   219   220   229   238   251   267 306];


accLeftSel = [125,85,86,73,106,80,81];
accRightSel = [357,286,287,255,295,281,283,282];




roiList = {'dACC', 'SPL', 'IPS', 'PFCl'}
% roiList = {'PFCl', 'PFClSal', 'PFClCtrl', 'dACC', 'accSal', 'accCtrl'}
subnetworkList = {'ContA', 'ContB', 'ContC', 'DorsAttnA', 'DorsAttnB', 'SalVenAttnA', 'SalVenAttnB', 'SomMotA', 'SalVenAttnB', 'VisualA', 'VisualB', 'VisualC'}
networkList = {'Cont', 'DorsAttn', 'SalVenAttn', 'SomMot', 'Visual'}



mkdir(fullfile(save_dir, 'regions'));
delete(fullfile(save_dir, 'regions', '*.nii'))



% t val
for nn = 1

    figure;
    tiledlayout('flow', 'padding', 'tight', 'TileSpacing','tight')
    c = load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/ScientificColourMaps7/vanimo/vanimo.mat');
    colormap(c.vanimo)

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
            case 'FEF'
                regSel = ismember(1:400, fefSel);
            case 'PFClSal'
                regSel = ismember(1:400, PFClSalSel);
            case 'PFClCtrl'
                regSel = ismember(1:400, PFClCtrlSel);
            case 'accSal'
                regSel = ismember(1:400, accSalSel);
            case 'accCtrl'
                regSel = ismember(1:400, accCtrlSel);
            case 'MOT'
                regSel = ismember(1:400, motSel);
            case 'VISA'
                regSel = ismember(1:400, visaSel);
            case 'M1Left'
                regSel = ismember(1:400, M1LeftSel);
            case 'dACCLeft'
                regSel = ismember(1:400, accLeftSel);
            case 'dACCRight'
                regSel = ismember(1:400, accRightSel);
            case 'PrCd'
                regSel = ismember(1:400, PrCdSel);
            case 'PFCd'
                regSel = ismember(1:400, PFCdSel);


            otherwise
                regSel = cellfun(@any, regexpi(regionTbl, regionList{rr}));
        end

        disp(regionList{rr})


        switch regionList{rr}
            case {'IPS', 'SPL'}

                clear r_tt r_td r_dt r_dd
                for pp = 1:npt

                    [r_tt(pp,1)] = corr(squeeze(R(1,1,regSel,pp)), squeeze(R(3,3,regSel,pp)), 'type', 'spearman');
                    [r_td(pp,1)] = corr(squeeze(R(1,1,regSel,pp)), squeeze(R(4,4,regSel,pp)), 'type', 'spearman');
                    [r_dt(pp,1)] = corr(squeeze(R(2,2,regSel,pp)), squeeze(R(3,3,regSel,pp)), 'type', 'spearman');
                    [r_dd(pp,1)] = corr(squeeze(R(2,2,regSel,pp)), squeeze(R(4,4,regSel,pp)), 'type', 'spearman');

                end

                %                 disp(regionList{rr})

                [~,pval,~, ~]=ttest([r_tt, r_td, r_dt, r_dd])

                %                 figure;
                %                 violinplot([r_tt, r_td, r_dt, r_dd]);


        end


        % print stats
        [h,pval,ci,stats] = ttest(reshape(squeeze(mean(R(:,:,regSel,:),3)), [nCond^2, npt])', [], 'alpha', .1);
        [hx,hy] = ind2sub([nCond, nCond], find(h));

        reP_t = reshape(stats.tstat, [nCond, nCond])
        reP_2tail = reshape(pval, [nCond, nCond])

        if ismember(regionList{rr}, {'SPL', 'IPS'})
            reP_d = reshape(stats.tstat ./ sqrt(stats.df+1), [nCond, nCond])
            reP_lower = reshape(ci(1,:), [nCond, nCond])
            reP_uppwr = reshape(ci(2,:), [nCond, nCond])
        end


        %         keyboard

        nexttile;
        imagesc(reshape(stats.tstat, [nCond, nCond]), [-5,5]);
        text(hx ,hy, '*', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', 12);

        title([regionList{rr}, ' T'])

        if rr == length(regionList)
            colorbar;
        end

        xticks(1:nCond); xticklabels(condLabel)
        yticks(1:nCond); yticklabels(condLabel)

        %         if rr == 1
        %         yticks(1:nCond); yticklabels(condLabel)
        %         end
        axis('square')






        % make mask
        parcelBrain = zeros(size(maskImg));
        parcelBrain(ismember(maskImg, find(regSel))) = 1;


        fname = fullfile(save_dir, 'regions', sprintf('%s.nii', regionList{rr}));

        % write
        VoOut      = struct(...
            'fname',    fname,...
            'dim',      maskVo.dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      maskVo.mat .* [-1 1 1 1]',...
            'n',        [1 1],...
            'descrip',  regionList{rr});

        spm_write_vol(VoOut, parcelBrain);




    end



end



%% directly compare SPL & IPS

clc;


accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];
fefSel = [315, 251, 116, 58];

PFClSalSel = [103, 104, 105, 307, 308, 309];
PFClCtrlSel = [59, 60, 61, 62, 252, 253, 254, 268,];
PFClCtrlSel = [278, 279, 280];

accSalSel = [    85    86   106   286   287   295];
accCtrlSel = [    73    80    81   255   281   282   283];


motSel = [176, 375];
visaSel =  [172   173   174   175   176   177   178   179   180   181   182   369   370   371   372   373   374   375   376   377   378   379   380];
M1LeftSel = [149, 150, 151];
PrCdSel = [   116   130   131   315   329   330   331];
PFCdSel = [     5    20    21    22    23    32    33    45    46    58    72    92   115   203   217   218   219   220   229   238   251   267 306];


accLeftSel = [125,85,86,73,106,80,81];
accRightSel = [357,286,287,255,295,281,283,282];



roiList = {'SPL', 'IPS'}






clear reg_resp reg_coh reg_cong


% t val
for nn = 1


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
            case 'FEF'
                regSel = ismember(1:400, fefSel);
            case 'PFClSal'
                regSel = ismember(1:400, PFClSalSel);
            case 'PFClCtrl'
                regSel = ismember(1:400, PFClCtrlSel);
            case 'accSal'
                regSel = ismember(1:400, accSalSel);
            case 'accCtrl'
                regSel = ismember(1:400, accCtrlSel);
            case 'MOT'
                regSel = ismember(1:400, motSel);
            case 'VISA'
                regSel = ismember(1:400, visaSel);
            case 'M1Left'
                regSel = ismember(1:400, M1LeftSel);
            case 'dACCLeft'
                regSel = ismember(1:400, accLeftSel);
            case 'dACCRight'
                regSel = ismember(1:400, accRightSel);
            case 'PrCd'
                regSel = ismember(1:400, PrCdSel);
            case 'PFCd'
                regSel = ismember(1:400, PFCdSel);


            otherwise
                regSel = cellfun(@any, regexpi(regionTbl, regionList{rr}));
        end


        reg_resp(:,:, rr) = [squeeze(mean(R(1,1,regSel,:),3)), squeeze(mean(R(2,2,regSel,:),3)), squeeze(mean(R(1,2,regSel,:),3))];
        reg_coh(:,:, rr) = [squeeze(mean(R(3,3,regSel,:),3)), squeeze(mean(R(4,4,regSel,:),3)), squeeze(mean(R(3,4,regSel,:),3))];
        reg_cong(:,:, rr) = [squeeze(mean(R(5,5,regSel,:),3))];


    end


    reg_list = roiList
    con_list = {'targ', 'dist', 'align'}

    for ii = 1:2

        for jj = 1:3

            % fprintf(['\n', reg_list{ii}, '\n'])
            % fprintf([con_list{jj}, '\n'])
            % fprintf('EVIDENCE\n')
            % 
            % 
            % [~,wn_p,wn_CI,wn_stats] = ttest(reg_resp(:,jj,ii))
            % wn_CI = wn_CI'
            % wn_bf = log10(bf_ttest('T', wn_stats.tstat, 'N', wn_stats.df+1))
            % 
            % 
            % fprintf(['\n', reg_list{ii}, '\n'])
            % fprintf([con_list{jj}, '\n'])
            % fprintf('COH\n')
            % 
            % [~,wn_p,wn_CI,wn_stats] = ttest(reg_coh(:,jj,ii))
            % wn_CI = wn_CI'
            % wn_bf = log10(bf_ttest('T', wn_stats.tstat, 'N', wn_stats.df+1))

        end


            fprintf(['\n', reg_list{ii}, '\n'])
            fprintf('CONG\n')

            [~,wn_p,wn_CI,wn_stats] = ttest(reg_cong(:,:,ii))
            wn_CI = wn_CI'
            wn_bf = log10(bf_ttest('T', wn_stats.tstat, 'N', wn_stats.df+1))


    end



    fprintf('CON RESP\n')

    [~,con_resp_pval,con_resp_CI,con_resp_stats] = ttest(reg_resp(:,:,1) - reg_resp(:,:,2))
    con_resp_CI = con_resp_CI'
    for ii = 1:3
        con_resp_bf = log10(bf_ttest('T', con_resp_stats.tstat(ii), 'N', con_resp_stats.df(ii)+1))
    end

    fprintf('CON COH\n')

    [~,con_coh_pval,con_coh_CI,con_coh_stats] = ttest(reg_coh(:,:,1) - reg_coh(:,:,2))
    con_coh_CI = con_coh_CI'
    for ii = 1:3
        con_coh_bf = log10(bf_ttest('T', con_coh_stats.tstat(ii), 'N', con_coh_stats.df(ii)+1))
    end


    fprintf('CON CONG\n')

    [~,con_cong_pval,con_cong_CI,con_cong_stats] = ttest(reg_cong(:,:,1) - reg_cong(:,:,2))
    con_cong_CI = con_cong_CI'
        con_cong_bf = log10(bf_ttest('T', con_cong_stats.tstat, 'N', con_cong_stats.df+1))
   


end














%% correlate betwen parcels within each region
accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];



% roiList = {'ExStrSup', 'SPL', 'IPS', 'IPL', 'PFCl', 'dACC'}
roiList = {'dACC','SPL', 'IPS'}
subnetworkList = {'ContA', 'ContB', 'ContC', 'DorsAttnA', 'DorsAttnB', 'SalVenAttnA', 'SalVenAttnB', 'SomMotA', 'SalVenAttnB', 'VisualA', 'VisualB', 'VisualC'}
networkList = {'Cont', 'DorsAttn', 'SalVenAttn', 'SomMot', 'Visual'}



mkdir(fullfile(save_dir, 'regions'));
delete(fullfile(save_dir, 'regions', '*.nii'))



% t val
for nn = 1

    %     figure;
    % tiledlayout('flow', 'padding', 'tight', 'TileSpacing','tight')

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




        % print stats

        disp(regionList{rr})

        r = [];
        for ii = 1:npt
            r(ii) = corr(squeeze(R(1,1,regSel,ii)), squeeze(R(2,2,regSel,ii)),'type', 'spearman');
        end
        %         nexttile; hold on;
        %         histogram(r)
        r_group = mean(r)
        [~,p_group] = ttest(r)

        [r_combined, p_combined] = corr(squeeze(mean(R(1,1,regSel,:),4)), squeeze(mean(R(2,2,regSel,:),4)), 'type', 'pearson')




    end



end


% == bootstrap

regionList = roiList;
regionTbl = f.info.region;

regSel1 = cellfun(@any, regexpi(regionTbl, regionList{2})); sum(regSel1)
regSel2 = cellfun(@any, regexpi(regionTbl, regionList{3})); sum(regSel2)

nsim = 1e5;
r_boot = nan(nsim,2);
for ii = 1:nsim

    bootPt = datasample(1:npt, npt);

    r_boot(ii,:) =  ...
        [(corr(squeeze(mean(R(1,1,regSel1,bootPt),4)), squeeze(mean(R(2,2,regSel1,bootPt),4)), 'type', 'pearson')), ...
        (corr(squeeze(mean(R(1,1,regSel2,bootPt),4)), squeeze(mean(R(2,2,regSel2,bootPt),4)), 'type', 'pearson'))];

end



figure;

% spl
disp('SPL')
[prctile(r_boot(:,1), 2.5), median(r_boot(:,1)), prctile(r_boot(:,1), 97.5)]

nexttile;
histogram(r_boot(:,1))
xline(prctile(r_boot(:,1), 2.5))
xline(prctile(r_boot(:,1), 97.5))
mean(r_boot(:,1)<0)
title('SPL')

%ips
disp('IPS')
[prctile(r_boot(:,2), 2.5), median(r_boot(:,2)), prctile(r_boot(:,2), 97.5)]
nexttile;
histogram(r_boot(:,2))
xline(prctile(r_boot(:,2), 2.5))
xline(prctile(r_boot(:,2), 97.5))
mean(r_boot(:,2)<0)
title('IPS')


%ips
disp('diff')
[(prctile(diff(r_boot,1,2), 2.5)), median(diff(r_boot,1,2)), (prctile(diff(r_boot,1,2), 97.5))]
nexttile;
histogram(diff(r_boot,1,2))
xline(prctile(diff(r_boot,1,2), 2.5))
xline(prctile(diff(r_boot,1,2), 97.5))
mean(diff(r_boot,1,2)>0)
title('diff')


% == bootstrap
%
% regionList = roiList;
% regionTbl = f.info.region;
%
% regSel1 = cellfun(@any, regexpi(regionTbl, regionList{2})); sum(regSel1)
% regSel2 = cellfun(@any, regexpi(regionTbl, regionList{3})); sum(regSel2)
%
%
% SPL = [squeeze(mean(R(1,1,regSel1,:),4)), squeeze(mean(R(2,2,regSel1,:),4))];
% IPS = [squeeze(mean(R(1,1,regSel2,:),4)), squeeze(mean(R(2,2,regSel2,:),4))];
%
% r_spl = corr(SPL(:,1), SPL(:,2))
% r_ips = corr(IPS(:,1), IPS(:,2))
%
%
%     r_boot(ii,:) =  ...
%         [atanh(corr(squeeze(mean(R(1,1,regSel1,bootPt),4)), squeeze(mean(R(2,2,regSel1,bootPt),4)), 'type', 'pearson')), ...
%         atanh(corr(squeeze(mean(R(1,1,regSel2,bootPt),4)), squeeze(mean(R(2,2,regSel2,bootPt),4)), 'type', 'pearson'))];
%
%
%
% r_boot = nan(10000,2);
% for ii = 1:10000
%
%     bootPt = datasample(1:npt, npt);
%
%     r_boot(ii,:) =  ...
%         [atanh(corr(squeeze(mean(R(1,1,regSel1,bootPt),4)), squeeze(mean(R(2,2,regSel1,bootPt),4)), 'type', 'pearson')), ...
%         atanh(corr(squeeze(mean(R(1,1,regSel2,bootPt),4)), squeeze(mean(R(2,2,regSel2,bootPt),4)), 'type', 'pearson'))];
%
% end







%% control network A vs C



accSel = [357,286,287,255,295,281,283,282,125,45,85,86,73,106,80,81];


roiList = {'ExStrSup', 'SPL', 'IPS', 'IPL', 'PFCl', 'dACC'}
subnetworkList = {'ContA', 'ContB', 'ContC', 'DorsAttnA', 'DorsAttnB', 'SalVenAttnA', 'SalVenAttnB', 'SomMotA', 'SalVenAttnB', 'VisualA', 'VisualB', 'VisualC'}
networkList = {'Cont', 'DorsAttn', 'SalVenAttn', 'SomMot', 'Visual'}



sqm = @(x,d) squeeze(mean(x,d));

if strcmp(analysisName, 'perf_parcel')


    figure;

    regionList = subnetworkList;
    regionTbl = f.info.subNetwork;


    regSel1 = cellfun(@any, regexpi(f.info.subNetwork, 'ContA'));
    regSel2 = cellfun(@any, regexpi(f.info.subNetwork, 'ContC'));




    % print stats
    rtCors = [...
        sqm(R(ismember(condLabel, 'rt'),ismember(condLabel, 'cohTarg'),regSel1,:),3),... ./ sqrt(sqm(R(5,5,regSel1,:),3).*sqm(R(3,3,regSel1,:),3)),...
        sqm(R(ismember(condLabel, 'rt'),ismember(condLabel, 'cohDist'),regSel1,:),3),... ./ sqrt(sqm(R(5,5,regSel1,:),3).*sqm(R(4,4,regSel1,:),3)),...
        sqm(R(ismember(condLabel, 'rt'),ismember(condLabel, 'cohTarg'),regSel2,:),3),... ./ sqrt(sqm(R(5,5,regSel2,:),3).*sqm(R(3,3,regSel2,:),3)),...
        sqm(R(ismember(condLabel, 'rt'),ismember(condLabel, 'cohDist'),regSel2,:),3)];% ./ sqrt(sqm(R(5,5,regSel2,:),3).*sqm(R(4,4,regSel2,:),3))];


    rtMeans = [...
        sqm(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'cohTarg'),regSel1,:),3),... ./ sqrt(sqm(R(5,5,regSel1,:),3).*sqm(R(3,3,regSel1,:),3)),...
        sqm(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'cohDist'),regSel1,:),3),... ./ sqrt(sqm(R(5,5,regSel1,:),3).*sqm(R(4,4,regSel1,:),3)),...
        sqm(R(ismember(condLabel, 'rt'),ismember(condLabel, 'rt'),regSel1,:),3),... ./ sqrt(sqm(R(5,5,regSel1,:),3).*sqm(R(4,4,regSel1,:),3)),...
        sqm(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'cohTarg'),regSel2,:),3),... ./ sqrt(sqm(R(5,5,regSel2,:),3).*sqm(R(3,3,regSel2,:),3)),...
        sqm(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'cohDist'),regSel2,:),3),...
        sqm(R(ismember(condLabel, 'rt'),ismember(condLabel, 'rt'),regSel2,:),3),...
        ];% ./ sqrt(sqm(R(5,5,regSel2,:),3).*sqm(R(4,4,regSel2,:),3))];



    rtdcors = [...
        (sqm(R(ismember(condLabel, 'rt'),ismember(condLabel, 'cohTarg'),regSel1,:),3)) ./ sqrt(mean(sqm(R(5,5,regSel1,:),3).*sqm(R(3,3,regSel1,:),3))),...
        (sqm(R(ismember(condLabel, 'rt'),ismember(condLabel, 'cohDist'),regSel1,:),3)) ./ sqrt(mean(sqm(R(5,5,regSel1,:),3).*sqm(R(4,4,regSel1,:),3))),...
        (sqm(R(ismember(condLabel, 'rt'),ismember(condLabel, 'cohTarg'),regSel2,:),3)) ./ sqrt(mean(sqm(R(5,5,regSel2,:),3).*sqm(R(3,3,regSel2,:),3))),...
        (sqm(R(ismember(condLabel, 'rt'),ismember(condLabel, 'cohDist'),regSel2,:),3)) ./ sqrt(mean(sqm(R(5,5,regSel2,:),3).*sqm(R(4,4,regSel2,:),3)))];







    [h,pval,~,stats] = ttest(rtCors(:,1)-rtCors(:,2), rtCors(:,3)-rtCors(:,4))



    nexttile; hold on;
    errorbar([1,2,3,4.5,5.5,6.5], nanmean(rtMeans), nanstd(rtMeans)./sqrt(sum(isfinite(rtMeans))), 'ok', 'MarkerFaceColor', 'auto', 'LineWidth',1.5)
    yline(0)
    xticks([1,2,3,4.5,5.5,6.5]); xticklabels({'ctrlA-targ','ctrlA-dist','ctrlA-RT','ctrlC-targ','ctrlC-dist','ctrlC-RT'})
    ylabel('r')
    title('feature encoding')



    nexttile; hold on;
    errorbar([1,2,3.5,4.5], nanmean(rtCors), nanstd(rtCors)./sqrt(sum(isfinite(rtCors))), 'ok', 'MarkerFaceColor', 'auto', 'LineWidth',1.5)
    xlim([.5,5])
    yline(0)
    xticks([1,2,3.5,4.5]); xticklabels({'ctrlA targ-rt','ctrlA dist-rt','ctrlC targ-rt','ctrlC dist-rt'})
    ylabel('r')
    title('salience - RT corr')



    nexttile; hold on;
    errorbar([1,2,3.5,4.5], nanmean(rtdcors), nanstd(rtdcors)./sqrt(sum(isfinite(rtdcors))), 'ok', 'MarkerFaceColor', 'auto', 'LineWidth',1.5)
    xlim([.5,5])
    yline(0)
    xticks([1,2,3.5,4.5]); xticklabels({'ctrlA targ-rt','ctrlA dist-rt','ctrlC targ-rt','ctrlC dist-rt'})
    ylabel('corrected r')
    title('salience - RT corr (noise corrected)')






    set(gca, 'TickDir', 'out', 'LineWidth', 1)
    ylabel('r')

end




%% behav correlations


if strcmp(analysisName, 'perf_parcel') | strcmp(analysisName, 'feature_parcel')


    accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];
    roiList = {'IPS', 'SPL', 'dACC'}


    d = load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/behavior/behav_table.mat');


    regionList = roiList;
    regionTbl = f.info.region;


    for rr = 1:length(regionList)

        switch regionList{rr}
            case 'dACC'
                regSel = ismember(1:400, accSel);
            otherwise
                regSel = cellfun(@any, regexpi(regionTbl, regionList{rr}));
        end

        cohTarg = squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'cohTarg'),regSel,:),3));
        cohDist = squeeze(mean(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'cohDist'),regSel,:),3));
        cohTargDist = squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'cohDist'),regSel,:),3));

        respTarg = squeeze(mean(R(ismember(condLabel, 'respTarg'),ismember(condLabel, 'respTarg'),regSel,:),3));
        respDist = squeeze(mean(R(ismember(condLabel, 'respDist'),ismember(condLabel, 'respDist'),regSel,:),3));
        respTargDist = squeeze(mean(R(ismember(condLabel, 'respTarg'),ismember(condLabel, 'respDist'),regSel,:),3));

        % add RSA to models

        t = d.tbl;
        [...
            t.cohTarg, t.cohDist, t.cohTargDist,...
            t.respTarg, t.respDist, t.respTargDist...
            ] = deal(t.pt*0);

        pts = unique(t.pt);
        center = @(x) x - nanmean(x);

        for pp = 1:length(pts)

            t.cohTarg(t.pt == pts(pp)) = cohTarg(pp);
            t.cohDist(t.pt == pts(pp)) = cohDist(pp);
            t.cohTargDist(t.pt == pts(pp)) = cohTargDist(pp);

            t.respTarg(t.pt == pts(pp)) = respTarg(pp);
            t.respDist(t.pt == pts(pp)) = respDist(pp);
            t.respTargDist(t.pt == pts(pp)) = respTargDist(pp);


        end

        t.cohTarg = center(t.cohTarg);
        t.cohDist = center(t.cohDist);
        t.cohTargDist2 =  center(abs(t.cohTargDist));

        t.respTarg = center(t.respTarg);
        t.respDist = center(t.respDist);
        t.respTargDist2 =  center(t.respTargDist);


        corr([...
            t.cohTarg, t.cohDist, t.cohTargDist,...
            t.respTarg, t.respDist, t.respTargDist...
            ])


        % fit models
        disp(regionList{rr})

        % coh RT
        rtMdl = fitlme(t, 'lrt ~ targ*cohTarg + dist*cohDist   + (1 + targ*cohTarg + dist*cohDist  | pt)',...
            'Exclude', ~(isfinite(t.rt) & t.acc == 1),...
            'FitMethod', 'REML', 'DummyVarCoding', 'effects');
        [~,~,ffx_rt] = fixedEffects(rtMdl, 'DFMethod', 'Satterthwaite')

        [p_coefrt, F_coefrt, r_coefrt] = coefTest(rtMdl, [0 0 0 0 0 1 1])



        % coh Acc
        accMdl = fitglme(t, 'acc ~ targ*cohTarg + dist*cohDist + (1 + targ*cohTarg + dist*cohDist  | pt)',...
            'Distribution', 'binomial', 'Exclude', ~(isfinite(t.rt)),...
            'FitMethod', 'Laplace', 'DummyVarCoding', 'effects');
        [~,~,ffx_acc] = fixedEffects(accMdl)

        [p_coefacc, F_coefacc, r_coefacc] = coefTest(accMdl, [0 0 0 0 0 1 1])



        %         % resp RT
        %         rtMdl = fitlme(t, 'lrt ~ targ + dist + targ:respTarg + dist:respDist + (targ + dist):respTargDist2  + (1 + targ + dist + targ:respTarg + dist:respDist + (targ + dist):respTargDist2  | pt)',...
        %             'Exclude', ~(isfinite(t.rt)),...
        %             'FitMethod', 'REML', 'DummyVarCoding', 'effects');
        %         [~,~,ffx_rt] = fixedEffects(rtMdl, 'DFMethod', 'Satterthwaite')


        %
        %
        %         % resp Acc
        %         accMdl = fitglme(t, 'acc ~ targ + dist + targ:respTarg + dist:respDist + (targ + dist):respTargDist2  + (1 + targ + dist + targ:respTarg + dist:respDist + (targ + dist):respTargDist2 | pt)',...
        %             'Distribution', 'binomial', 'Exclude', ~(isfinite(t.rt)),...
        %             'FitMethod', 'Laplace', 'DummyVarCoding', 'effects');
        %         [~,~,ffx_acc] = fixedEffects(accMdl)





    end



end




%% parcel-wide behavioral correlations


if strcmp(analysisName, 'perf_parcel') | strcmp(analysisName, 'feature_parcel')


    d = load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/behavior/behav_table.mat');

    cons = [find(ismember(condLabel, 'cohTarg')), find(ismember(condLabel, 'cohDist'))];
    fnames = {...
        '/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/parcelRSA_2022-06/perf_parcel/targRT.nii',...
        '/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/parcelRSA_2022-06/perf_parcel/distRT.nii'};


    [~,~,~,stats1] = ttest(squeeze(R(cons(1),cons(1),:,:))');
    [~,~,~,stats2] = ttest(squeeze(R(cons(2),cons(2),:,:))');

    regs = find(abs(stats1.tstat) > 3.4 & abs(stats2.tstat) > 3.4);
    [ffx_targ, ffx_dist] = deal([]);
    center = @(x) x - nanmean(x);


    parfor rr = 1:length(regs)

        fprintf('%d / %d\n', rr, length(regs)); % print

        regSel = regs(rr);

        cohTarg = squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'cohTarg'),regSel,:),3));
        cohDist = squeeze(mean(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'cohDist'),regSel,:),3));
        cohTargDist = squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'cohDist'),regSel,:),3));

        respTarg = squeeze(mean(R(ismember(condLabel, 'respTarg'),ismember(condLabel, 'respTarg'),regSel,:),3));
        respDist = squeeze(mean(R(ismember(condLabel, 'respDist'),ismember(condLabel, 'respDist'),regSel,:),3));
        respTargDist = squeeze(mean(R(ismember(condLabel, 'respTarg'),ismember(condLabel, 'respDist'),regSel,:),3));

        % add RSA to models

        t = d.tbl;
        [...
            t.cohTarg, t.cohDist, t.cohTargDist,...
            t.respTarg, t.respDist, t.respTargDist...
            ] = deal(t.pt*0);

        pts = unique(t.pt);

        for pp = 1:length(pts)

            t.cohTarg(t.pt == pts(pp)) = cohTarg(pp);
            t.cohDist(t.pt == pts(pp)) = cohDist(pp);
            t.cohTargDist(t.pt == pts(pp)) = cohTargDist(pp);

            t.respTarg(t.pt == pts(pp)) = respTarg(pp);
            t.respDist(t.pt == pts(pp)) = respDist(pp);
            t.respTargDist(t.pt == pts(pp)) = respTargDist(pp);


        end

        t.cohTarg = center(t.cohTarg);
        t.cohDist = center(t.cohDist);
        t.cohTargDist2 =  center(abs(t.cohTargDist));

        t.respTarg = center(t.respTarg);
        t.respDist = center(t.respDist);
        t.respTargDist2 =  center(t.respTargDist);


        % fit models
        %         disp(regionList{rr})

        % coh RT
        rtMdl = fitlme(t, 'lrt ~ targ*cohTarg + dist*cohDist + (1 + targ*cohTarg + dist*cohDist  | pt)',...
            'Exclude', ~(isfinite(t.rt)),...
            'FitMethod', 'REML', 'DummyVarCoding', 'effects');

        [~,~,ffx_model] = fixedEffects(rtMdl, 'DFMethod', 'Satterthwaite');
        ffx_targ(rr) = ffx_model.tStat(end-1);
        ffx_dist(rr) = ffx_model.tStat(end);

    end


    % reliablity R
    selBrain = maskImg;



    % print to parcellation
    parcelBrain = maskImg*0;
    for rr = 1:length(regs)
        parcelBrain(selBrain == regs(rr)) = ffx_targ(rr);
    end

    % write
    VoOut      = struct(...
        'fname',    fnames{1},...
        'dim',      maskVo.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      maskVo.mat .* [-1 1 1 1]',...
        'n',        [1 1],...
        'descrip',  'perf - targ');

    spm_write_vol(VoOut, parcelBrain);



    % print to parcellation
    parcelBrain = maskImg*0;
    for rr = 1:length(regs)
        parcelBrain(selBrain == regs(rr)) = ffx_dist(rr);
    end

    % write
    VoOut      = struct(...
        'fname',    fnames{2},...
        'dim',      maskVo.dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      maskVo.mat .* [-1 1 1 1]',...
        'n',        [1 1],...
        'descrip',  'perf - dist');

    spm_write_vol(VoOut, parcelBrain);




end






%% individual differences in feature-performance allignment


if strcmp(analysisName, 'perf_parcel')


    accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];
    fefSel = [315, 251, 116, 58];

    PFClSalSel = [103, 104, 105, 307, 308, 309];
    PFClCtrlSel = [59, 60, 61, 62, 252, 253, 254, 268, 278, 279, 280];

    accSalSel = [    85    86   106   286   287   295];
    accCtrlSel = [    73    80    81   255   281   282   283];


    motSel = [176, 375];
    visaSel =  [172   173   174   175   176   177   178   179   180   181   182   369   370   371   372   373   374   375   376   377   378   379   380];
    M1LeftSel = [149, 150, 151];
    PrCdSel = [   116   130   131   315   329   330   331];
    PFCdSel = [     5    20    21    22    23    32    33    45    46    58    72    92   115   203   217   218   219   220   229   238   251   267 306];


    accLeftSel = [125,85,86,73,106,80,81];
    accRightSel = [357,286,287,255,295,281,283,282];




    roiList = {'dACC', 'SPL', 'IPS', 'PFCl'}


    d = load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/behavior/targdist_rfx.mat');


    regionList = roiList;
    regionTbl = f.info.region;


    % figure;
    % tiledlayout(1,3, 'padding', 'tight', 'TileSpacing','tight')

    for rr = 1:length(regionList)

        switch regionList{rr}
            case 'dACC'
                regSel = ismember(1:400, accSel);
            case 'FEF'
                regSel = ismember(1:400, fefSel);
            case 'PFClSal'
                regSel = ismember(1:400, PFClSalSel);
            case 'PFClCtrl'
                regSel = ismember(1:400, PFClCtrlSel);
            case 'accSal'
                regSel = ismember(1:400, accSalSel);
            case 'accCtrl'
                regSel = ismember(1:400, accCtrlSel);
            case 'MOT'
                regSel = ismember(1:400, motSel);
            case 'VISA'
                regSel = ismember(1:400, visaSel);
            case 'M1Left'
                regSel = ismember(1:400, M1LeftSel);
            case 'dACCLeft'
                regSel = ismember(1:400, accLeftSel);
            case 'dACCRight'
                regSel = ismember(1:400, accRightSel);
            case 'PrCd'
                regSel = ismember(1:400, PrCdSel);
            case 'PFCd'
                regSel = ismember(1:400, PFCdSel);


            otherwise
                regSel = cellfun(@any, regexpi(regionTbl, regionList{rr}));
        end


        disp(regionList{rr})
        % == corr
        [r, p] = corr([...
            squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'rt'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'acc'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'rt'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'acc'),regSel,:),3)),...
            ]);

        p

        con = (r(1,4) + r(2,3)) - (r(1,3) + r(2,4))



        [hx,hy] = ind2sub(size(r), find(p<.05));


        % == partial cor

        % targ
        ctrlmx = [squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'cohTarg'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'rt'),ismember(condLabel, 'rt'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'acc'),ismember(condLabel, 'acc'),regSel,:),3)),...
            ];

        cormx = [...
            squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'rt'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'acc'),regSel,:),3)),...
            ];


        [targ_partial_r, targ_partial_p] = partialcorr(cormx, ctrlmx)


        % dist
        ctrlmx = [squeeze(mean(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'cohDist'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'rt'),ismember(condLabel, 'rt'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'acc'),ismember(condLabel, 'acc'),regSel,:),3)),...
            ];


        cormx = [...
            squeeze(mean(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'rt'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'acc'),regSel,:),3)),...
            ];


        [dist_partial_r, dist_partial_p] = partialcorr(cormx, ctrlmx)




        % == corr
        [~,enc_p,enc_ci,enc_stats] = ttest([...
            squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'rt'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'acc'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'rt'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'acc'),regSel,:),3)),...
            ]);

        [enc_hx,enc_hy] = ind2sub(size(r), find(diag(enc_p)<.1 & diag(enc_p)>0));

        enc_p
        enc_ci = enc_ci'






        %         nexttile;

        %         imagesc(r, [-1, 1]);
        %         title(regionList{rr})
        %         xticks(1:4); xticklabels({'targRT', 'targAcc', 'distRT', 'distAcc'})
        %         yticks(1:4); yticklabels({'targRT', 'targAcc', 'distRT', 'distAcc'})
        %         text(hx ,hy, '*', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', 12);
        %         text(enc_hx ,enc_hy, '*', 'color', 'b', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', 12);
        %
        %         axis('square')


        fg=figure;
        corrplot([...
            squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'rt'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'cohTarg'),ismember(condLabel, 'acc'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'rt'),regSel,:),3)),...
            squeeze(mean(R(ismember(condLabel, 'cohDist'),ismember(condLabel, 'acc'),regSel,:),3)),...
            ],  'TestR', 'on', 'Type', 'Spearman')
        title(regionList{rr})
        fg.Units = 'inches';
        fg.OuterPosition = [0 0 6 6.5];



    end



end






%% performance encoding across regions


% if strcmp(analysisName, 'perf_parcel')


accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];
roiList = {'dACC', 'SPL', 'IPS'}


d = load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/behavior/targdist_rfx.mat');


regionList = roiList;
regionTbl = f.info.region;



clear regRt regAcc
for rr = 1:length(regionList)

    switch regionList{rr}
        case 'dACC'
            regSel = ismember(1:400, accSel);
        otherwise
            regSel = cellfun(@any, regexpi(regionTbl, regionList{rr}));
    end



    % == corr
    regRt(:,rr) = squeeze(mean(R(ismember(condLabel, 'cong'),ismember(condLabel, 'cong'),regSel,:),3));
    regAcc(:,rr) = squeeze(mean(R(ismember(condLabel, 'cong'),ismember(condLabel, 'cong'),regSel,:),3));


end

regRt = regRt - mean(regRt,2);
regAcc = regAcc - mean(regAcc,2);


figure;
nexttile;hold on;
errorbar(mean(regRt), std(regRt)./sqrt(npt), 'or', 'MarkerFaceColor','r')
xlim([.5, 3.5])

nexttile;hold on;
errorbar(mean(regAcc), std(regAcc)./sqrt(npt), 'ob', 'MarkerFaceColor','b')
xlim([.5, 3.5])






% end






%% MDS


addpath(genpath('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/pcm_toolbox'));
accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];

load('bam.mat')
cmap = bam;
cidx = round(linspace(1, 256, 8));


roiList = {'dACC', 'SPL', 'IPS', 'PFCl'}
subnetworkList = {'ContA', 'ContB', 'ContC', 'DorsAttnA', 'DorsAttnB', 'SalVenAttnA', 'SalVenAttnB', 'SomMotA', 'SalVenAttnB', 'VisualA', 'VisualB', 'VisualC'}
networkList = {'Cont', 'DorsAttn', 'SalVenAttn', 'SomMot', 'Visual'}


if strcmp(analysisName, 'margCohTask_parcel')

    % t val
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




            % print stats
            M = mean(reshape(squeeze(mean(G(:,:,regSel,:),3)), [nCond^2, npt])');
            rM = reshape(M, [nCond, nCond]);

            %                         rM = sqrt(diag([ones(1,8), ones(1,4)*2]/1.3333))*rM*sqrt(diag([ones(1,8), ones(1,4)*2]/1.3333));
            [L,l] = pcm_classicalMDS(rM([1:4,9:12], [1:4,9:12]));



            % combined
            mds = L;
            nexttile; hold on;

            plot(mds(1:4,1), mds(1:4,2),  'or', 'MarkerSize', 8); lsline;
            plot(mds(5:8,1), mds(5:8,2),  'og', 'MarkerSize', 8);


            nexttile; hold on;

            plot(mds(1:4,2), mds(1:4,3),  'or', 'MarkerSize', 8); lsline;
            plot(mds(5:8,2), mds(5:8,3),  'og', 'MarkerSize', 8);

            %             for ii = 1:4
            %                 plot(mds(ii,1), mds(ii,2),  'o',...
            %                     'Color', cmap(cidx(ii),:), 'MarkerFaceColor', cmap(cidx(ii),:), 'MarkerSize', 10, 'MarkerEdgeColor','k');
            %             end
            %
            %             for ii = 1:4
            %                 plot(mds(8+ii,1), mds(8+ii,2), 'd',...
            %                     'Color', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerFaceColor', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerSize', 10, 'MarkerEdgeColor','k');
            %             end


            title([regionList{rr}])
            set(gca, 'TickDir', 'out', 'LineWidth', 1)




        end





        %         for rr = 1:length(regionList)
        %
        %             switch regionList{rr}
        %                 case 'dACC'
        %                     regSel = ismember(1:400, accSel);
        %                 otherwise
        %                     regSel = cellfun(@any, regexpi(regionTbl, regionList{rr}));
        %             end
        %
        %
        %
        %
        %             % print stats
        %             M = mean(reshape(squeeze(mean(G(:,:,regSel,:),3)), [nCond^2, npt])');
        %             rM = reshape(M, [nCond, nCond]);
        %
        %             %             rM = sqrt(diag([ones(1,8), ones(1,4)*2]/1.3333))*rM*sqrt(diag([ones(1,8), ones(1,4)*2]/1.3333));
        %             [L,l] =pcm_classicalMDS(rM(1:8, 1:8));
        %
        %
        %
        %             % combined
        %             mds = L;
        %             nexttile; hold on;
        %
        %             plot(mds(1:4,1), mds(1:4,2),  'or', 'MarkerSize', 1);
        %             plot(mds(5:8,1), mds(5:8,2),  'og', 'MarkerSize', 1);
        %
        %             for ii = 1:8
        %                 plot(mds(ii,1), mds(ii,2),  'o',...
        %                     'Color', cmap(cidx(ii),:), 'MarkerFaceColor', cmap(cidx(ii),:), 'MarkerSize', 10, 'MarkerEdgeColor','k');
        %             end
        %
        %             %             for ii = 1:4
        %             %                 plot(mds(8+ii,1), mds(8+ii,2), 'd',...
        %             %                     'Color', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerFaceColor', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerSize', 10);
        %             %             end
        %
        %
        %             title([regionList{rr}, ' target-only'])
        %             set(gca, 'TickDir', 'out', 'LineWidth', 1)
        %
        %
        %
        %
        %         end








    end


end







%% MDS WITH RESP


addpath(genpath('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/pcm_toolbox'));
accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];
fefSel = [315, 251, 116, 58];
PFClSalSel = [103, 104, 105, 307, 308, 309];
PFClCtrlSel = [59, 60, 61, 62, 252, 253, 254, 268, 278, 279, 280];
accSalSel = [    85    86   106   286   287   295];
accCtrlSel = [    73    80    81   255   281   282   283];
motSel = [176, 375];
visaSel =  [172   173   174   175   176   177   178   179   180   181   182   369   370   371   372   373   374   375   376   377   378   379   380];
M1LeftSel = [149, 150, 151];
PrCdSel = [   116   130   131   315   329   330   331];
PFCdSel = [     5    20    21    22    23    32    33    45    46    58    72    92   115   203   217   218   219   220   229   238   251   267 306];


accLeftSel = [125,85,86,73,106,80,81];
accRightSel = [357,286,287,255,295,281,283,282];



load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/ScientificColourMaps7/roma/roma.mat')
cmap = roma;
cidx = round(linspace(1, 256, 8));


roiList = {'dACC', 'PFCl', 'SPL', 'IPS'}
% roiList = {'dACC', 'PFCl', 'SPL', 'IPS'}
subnetworkList = {'ContA', 'ContB', 'ContC', 'DorsAttnA', 'DorsAttnB', 'SalVenAttnA', 'SalVenAttnB', 'SomMotA', 'SalVenAttnB', 'VisualA', 'VisualB', 'VisualC'}
networkList = {'Cont', 'DorsAttn', 'SalVenAttn', 'SomMot', 'Visual'}


if strcmp(analysisName, 'margResp_parcel')

    % t val
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
                case 'FEF'
                    regSel = ismember(1:400, fefSel);
                case 'PFClSal'
                    regSel = ismember(1:400, PFClSalSel);
                case 'PFClCtrl'
                    regSel = ismember(1:400, PFClCtrlSel);
                case 'accSal'
                    regSel = ismember(1:400, accSalSel);
                case 'accCtrl'
                    regSel = ismember(1:400, accCtrlSel);
                case 'MOT'
                    regSel = ismember(1:400, motSel);
                case 'VISA'
                    regSel = ismember(1:400, visaSel);
                case 'M1Left'
                    regSel = ismember(1:400, M1LeftSel);
                case 'dACCLeft'
                    regSel = ismember(1:400, accLeftSel);
                case 'dACCRight'
                    regSel = ismember(1:400, accRightSel);
                case 'PrCd'
                    regSel = ismember(1:400, PrCdSel);
                case 'PFCd'
                    regSel = ismember(1:400, PFCdSel);
                otherwise
                    regSel = cellfun(@any, regexpi(regionTbl, regionList{rr}));
            end




            % print stats
            M = mean(reshape(squeeze(mean(R(:,:,regSel,:),3)), [nCond^2, npt])');
            rM = reshape(M, [nCond, nCond]);

            %                         rM = sqrt(diag([ones(1,8), ones(1,4)*2]/1.3333))*rM*sqrt(diag([ones(1,8), ones(1,4)*2]/1.3333));


            [L,l] = pcm_classicalMDS(rM);


            % combined
            mds = L;
            nexttile; hold on;

            plot(mds(1:4,1), mds(1:4,2),  'or', 'MarkerSize', 1); lsline;
            plot(mds(5:8,1), mds(5:8,2),  'og', 'MarkerSize', 1); lsline;

            for ii = 1:8
                plot(mds(ii,1), mds(ii,2),  'o',...
                    'Color', cmap(cidx(ii),:), 'MarkerFaceColor', cmap(cidx(ii),:), 'MarkerSize', 10, 'MarkerEdgeColor','k');
            end

            for ii = 1:4
                plot(mds(8+ii,1), mds(8+ii,2), 'd',...
                    'Color', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerFaceColor', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerSize', 10, 'MarkerEdgeColor','k');
            end



            title([regionList{rr}])
            set(gca, 'TickDir', 'out', 'LineWidth', 1)
            axis('square')




        end


    end









    % higher dim
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




            % print stats
            M = mean(reshape(squeeze(mean(R(:,:,regSel,:),3)), [nCond^2, npt])');
            rM = reshape(M, [nCond, nCond]);

            %                         rM = sqrt(diag([ones(1,8), ones(1,4)*2]/1.3333))*rM*sqrt(diag([ones(1,8), ones(1,4)*2]/1.3333));
            [L,l] = pcm_classicalMDS(rM);



            % combined
            mds = L;
            nexttile; hold on;

            %             plot(mds(1:4,2), mds(1:4,3),  'or', 'MarkerSize', 1); lsline;
            %             plot(mds(5:8,2), mds(5:8,3),  'og', 'MarkerSize', 1); lsline;
            %
            %
            %             for ii = 1:8
            %                 plot(mds(ii,2), mds(ii,3),  'o',...
            %                     'Color', cmap(cidx(ii),:), 'MarkerFaceColor', cmap(cidx(ii),:), 'MarkerSize', 10, 'MarkerEdgeColor','k');
            %             end
            %
            %             for ii = 1:4
            %                 plot(mds(8+ii,2), mds(8+ii,3), 'd',...
            %                     'Color', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerFaceColor', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerSize', 10, 'MarkerEdgeColor','k');
            %             end

            %
            %              plot3(mds(1:4,1), mds(1:4,2), mds(1:4,3),  'or', 'MarkerSize', 1); lsline;
            %              plot3(mds(5:8,1), mds(5:8,2), mds(5:8,3),  'or', 'MarkerSize', 1); lsline;
            %             plot3(mds(10:12,1), mds(10:12,2), mds(10:12,3),  'og', 'MarkerSize', 1); lsline;
            %
            %
            % plot3(...
            %     [mean(mds([4,8],1)), mean(mds([1,5],1))]',...
            %     [mean(mds([4,8],2)), mean(mds([1,5],2))]',...
            %     [mean(mds([4,8],3)), mean(mds([1,5],3))]'...
            %     )
            %
            %





            % plot lines

            plot3(...
                mds(1:4,1),...
                mds(1:4,2),...
                mds(1:4,3), '-', 'color', cmap(cidx(1),:),...
                'LineWidth', 2 ...
                )

            plot3(...
                mds(5:8,1),...
                mds(5:8,2),...
                mds(5:8,3), '-', 'color', cmap(cidx(end),:),...
                'LineWidth', 2 ...
                )


            plot3(...
                [mean(mds([9,12],1)), mean(mds([10,11],1))],...
                [mean(mds([9,12],2)), mean(mds([10,11],2))],...
                [mean(mds([9,12],3)), mean(mds([10,11],3))],...
                '-k', 'LineWidth', 2 ...
                )





            for ii = 1:8
                plot3(mds(ii,1), mds(ii,2), mds(ii,3),  'o',...
                    'Color', cmap(cidx(ii),:), 'MarkerFaceColor', cmap(cidx(ii),:), 'MarkerSize', 12, 'MarkerEdgeColor','k');

                %                 plot3(mds(ii,1), zeros(size(mds(ii,1))) + .3, mds(ii,3),'o', 'MarkerSize', 10,'color',[0.8 0.8 0.8]);

            end

            for ii = 1:4
                plot3(mds(8+ii,1), mds(8+ii,2), mds(8+ii,3), 'd',...
                    'Color', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerFaceColor', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerSize', 12, 'MarkerEdgeColor','k');

                %                 plot3(mds(8+ii,1), zeros(size(mds(8+ii,1))) + .3, mds(ii,3),'o', 'MarkerSize', 10,'color',[0.8 0.8 0.8]);

            end








            % shadow
            xl = xlim*1.1;
            yl = ylim*1.1;
            zl = zlim*1.1;

            plot3(...
                mds(1:8,1),...
                mds(1:8,2),...
                0*mds(1:8,3) + zl(1), 'o', ...
                'color', [.8,.8,.8],...
                'MarkerSize', 14, ...
                'MarkerFaceColor', [.8,.8,.8]...
                )


             plot3(...
                mds(9:12,1),...
                mds(9:12,2),...
                0*mds(9:12,3) + zl(1), 'd', ...
                'color', [.8,.8,.8],...
                'MarkerSize', 14, ...
                'MarkerFaceColor', [.8,.8,.8]...
                )


            % z project
            plot3(...
                mds(1:4,1),...
                mds(1:4,2),...
                0*mds(1:4,3) + zl(1), '-', ...
                'LineWidth', 2.5,...
                'color', [.8,.8,.8],...
                'MarkerFaceColor', [.8,.8,.8]...
                )

            plot3(...
                mds(5:8,1),...
                mds(5:8,2),...
                0*mds(5:8,3) + zl(1), '-', ...
                'LineWidth', 2.5,...
                'color', [.8,.8,.8],...
                'MarkerFaceColor', [.8,.8,.8]...
                )


            plot3(...
                [mean(mds([9,12],1)), mean(mds([10,11],1))],...
                [mean(mds([9,12],2)), mean(mds([10,11],2))],...
                [0*mean(mds([9,12],3)) + zl(1), 0*mean(mds([10,11],3)) + zl(1)], '-', ...
                'LineWidth', 2.5,...
                'color', [.8,.8,.8],...
                'MarkerFaceColor', [.8,.8,.8]...
                )



            % y project
            plot3(...
                mds(1:4,1),...
                0*mds(1:4,2) + yl(2),...
                mds(1:4,3), '-', ...
                'LineWidth', 2,...
                'color', [.8,.8,.8],...
                'MarkerFaceColor', [.8,.8,.8]...
                )

            plot3(...
                mds(5:8,1),...
                0*mds(5:8,2) + yl(2),...
                mds(5:8,3), '-', ...
                'LineWidth', 2,...
                'color', [.8,.8,.8],...
                'MarkerFaceColor', [.8,.8,.8]...
                )


            plot3(...
                [mean(mds([9,12],1)), mean(mds([10,11],1))],...
                [0*mean(mds([9,12],2)) + yl(2), 0*mean(mds([10,11],2))+ yl(2)],...
                [mean(mds([9,12],3)), mean(mds([10,11],3))], '-', ...
                'LineWidth', 2,...
                'color', [.8,.8,.8],...
                'MarkerFaceColor', [.8,.8,.8]...
                )





            % x project
            plot3(...
                0*mds(1:4,1) + xl(2),...
                mds(1:4,2),...
                mds(1:4,3), '-', ...
                'LineWidth', 2,...
                'color', [.8,.8,.8],...
                'MarkerFaceColor', [.8,.8,.8]...
                )

            plot3(...
                0*mds(5:8,1) + xl(2),...
                mds(5:8,2),...
                mds(5:8,3), '-', ...
                'LineWidth', 2,...
                'color', [.8,.8,.8],...
                'MarkerFaceColor', [.8,.8,.8]...
                )


            plot3(...
                [0*mean(mds([9,12],1)) + xl(2), 0*mean(mds([10,11],1))+ xl(2)],...
                [mean(mds([9,12],2)), mean(mds([10,11],2))],...
                [mean(mds([9,12],3)), mean(mds([10,11],3))], '-', ...
                'LineWidth', 2,...
                'color', [.8,.8,.8],...
                'MarkerFaceColor', [.8,.8,.8]...
                )

         
    
    
    %             title([regionList{rr}])
                set(gca, 'TickDir', 'None', 'LineWidth', 1.5)
                axis('equal')
                grid on
    
                xlim(xl)
                ylim(yl)
                zlim(zl)
                set(gca, 'CameraPosition', [-1.1558,   -1.6413,    1.0410])
                xlabel('Eigenvector 1', 'FontName', 'Myriad Pro', 'FontSize',14)
                ylabel('Eigenvector 2', 'FontName', 'Myriad Pro', 'FontSize',14)
                zlabel('Eigenvector 3', 'FontName', 'Myriad Pro', 'FontSize',14)
                xticklabels([])
                yticklabels([])
                zticklabels([])

            %             ylim([-.1, .1])


        end


    end



% 
% 
%     % t val
%     for nn = 1
% 
%         figure;
%         tiledlayout('flow', 'padding', 'compact', 'TileSpacing','compact')
% 
% 
%         switch nn
% 
%             case 1
% 
%                 regionList = roiList;
%                 regionTbl = f.info.region;
% 
%             case 2
% 
%                 regionList = subnetworkList;
%                 regionTbl = f.info.subNetwork;
% 
%             case 3
% 
%                 regionList = networkList;
%                 regionTbl = f.info.network;
% 
%         end
% 
% 
%         for rr = 1:length(regionList)
% 
%             switch regionList{rr}
%                 case 'dACC'
%                     regSel = ismember(1:400, accSel);
%                 case 'FEF'
%                     regSel = ismember(1:400, fefSel);
%                 case 'PFClSal'
%                     regSel = ismember(1:400, PFClSalSel);
%                 case 'PFClCtrl'
%                     regSel = ismember(1:400, PFClCtrlSel);
%                 case 'accSal'
%                     regSel = ismember(1:400, accSalSel);
%                 case 'accCtrl'
%                     regSel = ismember(1:400, accCtrlSel);
%                 case 'MOT'
%                     regSel = ismember(1:400, motSel);
%                 case 'VISA'
%                     regSel = ismember(1:400, visaSel);
%                 case 'M1Left'
%                     regSel = ismember(1:400, M1LeftSel);
%                 case 'dACCLeft'
%                     regSel = ismember(1:400, accLeftSel);
%                 case 'dACCRight'
%                     regSel = ismember(1:400, accRightSel);
%                 case 'PrCd'
%                     regSel = ismember(1:400, PrCdSel);
%                 case 'PFCd'
%                     regSel = ismember(1:400, PFCdSel);
%                 otherwise
%                     regSel = cellfun(@any, regexpi(regionTbl, regionList{rr}));
%             end
% 
% 
% 
% 
%             % print stats
%             M = mean(reshape(squeeze(mean(R(:,:,regSel,:),3)), [nCond^2, npt])');
%             rM = reshape(M, [nCond, nCond]);
% 
%             %                         rM = sqrt(diag([ones(1,8), ones(1,4)*2]/1.3333))*rM*sqrt(diag([ones(1,8), ones(1,4)*2]/1.3333));
% 
% 
%             [L,l] = pcm_classicalMDS(rM);
% 
% 
%             % combined
%             nexttile; hold on;
% 
%             plot(cumsum(l(l>0))/sum(l(l>0)), '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'k')
% 
% 
%             title([regionList{rr}])
%             set(gca, 'TickDir', 'out', 'LineWidth', 1)
%             axis('square')
% 
% 
% 
% 
%         end
% 
% 
%     end











end










%% MDS WITH RESP (3D)


addpath(genpath('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/pcm_toolbox'));
accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];
fefSel = [315, 251, 116, 58];
PFClSalSel = [103, 104, 105, 307, 308, 309];
PFClCtrlSel = [59, 60, 61, 62, 252, 253, 254, 268, 278, 279, 280];
accSalSel = [    85    86   106   286   287   295];
accCtrlSel = [    73    80    81   255   281   282   283];
motSel = [176, 375];
visaSel =  [172   173   174   175   176   177   178   179   180   181   182   369   370   371   372   373   374   375   376   377   378   379   380];
M1LeftSel = [149, 150, 151];
PrCdSel = [   116   130   131   315   329   330   331];
PFCdSel = [     5    20    21    22    23    32    33    45    46    58    72    92   115   203   217   218   219   220   229   238   251   267 306];


accLeftSel = [125,85,86,73,106,80,81];
accRightSel = [357,286,287,255,295,281,283,282];



load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/ScientificColourMaps7/roma/roma.mat')
cmap = roma;
cidx = round(linspace(1, 256, 8));


roiList = {'dACC', 'PFCl', 'SPL', 'IPS'}
% roiList = {'dACC', 'PFCl', 'SPL', 'IPS'}
subnetworkList = {'ContA', 'ContB', 'ContC', 'DorsAttnA', 'DorsAttnB', 'SalVenAttnA', 'SalVenAttnB', 'SomMotA', 'SalVenAttnB', 'VisualA', 'VisualB', 'VisualC'}
networkList = {'Cont', 'DorsAttn', 'SalVenAttn', 'SomMot', 'Visual'}


if strcmp(analysisName, 'margResp_parcel')





        fg=figure;
 
        tiledlayout('flow', 'padding', 'compact', 'TileSpacing','compact')




    % higher dim
    for nn = 1

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




            % print stats
            M = mean(reshape(squeeze(mean(R(:,:,regSel,:),3)), [nCond^2, npt])');
            rM = reshape(M, [nCond, nCond]);




            [L,l] = pcm_classicalMDS(rM);

                                    nexttile; imagesc(rM, [-.02, .02]); title(roiList{rr}); colormap(roma)
                        nexttile; imagesc(L, [-.1, .1]); title(roiList{rr}); colormap(roma)



  
%             [L,l] = pcm_classicalMDS(rM,...
%                 'contrast', [...
%                 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0;...
%                 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1;...
%                 0.5, 0.1667, -0.1667, -0.5000, -0.5000, -0.1667, 0.1667, 0.5000, 0,0,0,0]');



            % combined



            mds = L;
            nexttile; hold on;


            % plot lines

            plot3(...
                mds(1:4,1),...
                mds(1:4,2),...
                mds(1:4,3), '-', 'color', cmap(cidx(1),:),...
                'LineWidth', 3 ...
                )

            plot3(...
                mds(5:8,1),...
                mds(5:8,2),...
                mds(5:8,3), '-', 'color', cmap(cidx(end),:),...
                'LineWidth', 3 ...
                )


            plot3(...
                [mean(mds([9,12],1)), mean(mds([10,11],1))],...
                [mean(mds([9,12],2)), mean(mds([10,11],2))],...
                [mean(mds([9,12],3)), mean(mds([10,11],3))],...
                '-k', 'LineWidth', 3 ...
                )


            for ii = 1:8
                plot3(mds(ii,1), mds(ii,2), mds(ii,3),  'o',...
                    'Color', cmap(cidx(ii),:), 'MarkerFaceColor', cmap(cidx(ii),:), 'MarkerSize', 18, 'MarkerEdgeColor','k');

                %                 plot3(mds(ii,1), zeros(size(mds(ii,1))) + .3, mds(ii,3),'o', 'MarkerSize', 10,'color',[0.8 0.8 0.8]);

            end

            for ii = 1:4
                plot3(mds(8+ii,1), mds(8+ii,2), mds(8+ii,3), 'd',...
                    'Color', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerFaceColor', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerSize', 18, 'MarkerEdgeColor','k');

                %                 plot3(mds(8+ii,1), zeros(size(mds(8+ii,1))) + .3, mds(ii,3),'o', 'MarkerSize', 10,'color',[0.8 0.8 0.8]);

            end








            % shadow
            xl = xlim*1.1;
            yl = ylim*1.1;
            zl = zlim*1.1;

            shadowCol = [.75, .75, .75];

            zcol1 =  0*cmap(cidx(1),:) + .85*[.85, .85, .85];
            zcol2 =  0*cmap(cidx(end),:) + .85*[.85, .85, .85];
            zcol3 = 0*[0,0,0] + .85*[.85, .85, .85];


            col1 =  .15*cmap(cidx(1),:) + .85*[.9, .9, .9];
            col2 =  .15*cmap(cidx(end),:) + .85*[.9, .9, .9];
            col3 = .15*[0,0,0] + .85*[.9, .9, .9];


            plot3(...
                mds(1:4,1),...
                mds(1:4,2),...
                0*mds(1:4,3) + zl(1), 'o', ...
                'color', zcol1,...
                'MarkerSize', 18, ...
                'MarkerFaceColor', zcol1...
                )


            plot3(...
                mds(5:8,1),...
                mds(5:8,2),...
                0*mds(5:8,3) + zl(1), 'o', ...
                'color', zcol2,...
                'MarkerSize', 18, ...
                'MarkerFaceColor', zcol2...
                )


             plot3(...
                mds(9:10,1),...
                mds(9:10,2),...
                0*mds(9:10,3) + zl(1), 'd', ...
                'color', zcol1,...
                'MarkerSize', 18, ...
                'MarkerFaceColor', zcol1...
                )

             plot3(...
                 mds(11:12,1),...
                 mds(11:12,2),...
                 0*mds(11:12,3) + zl(1), 'd', ...
                 'color', zcol2,...
                 'MarkerSize', 18, ...
                 'MarkerFaceColor', zcol2...
                 )


            % z project
            plot3(...
                mds(1:4,1),...
                mds(1:4,2),...
                0*mds(1:4,3) + zl(1), '-', ...
                'LineWidth', 3,...
                'color', zcol1,...
                'MarkerFaceColor', zcol1...
                )

            plot3(...
                mds(5:8,1),...
                mds(5:8,2),...
                0*mds(5:8,3) + zl(1), '-', ...
                'LineWidth', 3,...
                'color', zcol2,...
                'MarkerFaceColor', zcol2...
                )


            plot3(...
                [mean(mds([9,12],1)), mean(mds([10,11],1))],...
                [mean(mds([9,12],2)), mean(mds([10,11],2))],...
                [0*mean(mds([9,12],3)) + zl(1), 0*mean(mds([10,11],3)) + zl(1)], '-', ...
                'LineWidth', 3,...
                'color', zcol3,...
                'MarkerFaceColor', zcol3...
                )



            % y project
            plot3(...
                mds(1:4,1),...
                0*mds(1:4,2) + yl(2),...
                mds(1:4,3), '-', ...
                'LineWidth', 2,...
                'color', col1,...
                'MarkerFaceColor', col1...
                )

            plot3(...
                mds(5:8,1),...
                0*mds(5:8,2) + yl(2),...
                mds(5:8,3), '-', ...
                'LineWidth', 2,...
                'color', col2,...
                'MarkerFaceColor', col2...
                )


            plot3(...
                [mean(mds([9,12],1)), mean(mds([10,11],1))],...
                [0*mean(mds([9,12],2)) + yl(2), 0*mean(mds([10,11],2))+ yl(2)],...
                [mean(mds([9,12],3)), mean(mds([10,11],3))], '-', ...
                'LineWidth', 2,...
                'color', col3,...
                'MarkerFaceColor', col3...
                )





            % x project
            plot3(...
                0*mds(1:4,1) + xl(2),...
                mds(1:4,2),...
                mds(1:4,3), '-', ...
                'LineWidth', 2,...
                'color', col1,...
                'MarkerFaceColor', col1...
                )

            plot3(...
                0*mds(5:8,1) + xl(2),...
                mds(5:8,2),...
                mds(5:8,3), '-', ...
                'LineWidth', 2,...
                'color', col2,...
                'MarkerFaceColor', col2...
                )


            plot3(...
                [0*mean(mds([9,12],1)) + xl(2), 0*mean(mds([10,11],1))+ xl(2)],...
                [mean(mds([9,12],2)), mean(mds([10,11],2))],...
                [mean(mds([9,12],3)), mean(mds([10,11],3))], '-', ...
                'LineWidth', 2,...
                'color', col3,...
                'MarkerFaceColor', col3...
                )




            %             title([regionList{rr}])
            set(gca, 'TickDir', 'None', 'LineWidth', 1.5)
            axis('equal')
            grid on;

            xlim(xl)
            ylim(yl)
            zlim(zl)
            xlabel('Eigenvector 1', 'FontName', 'Myriad Pro', 'FontSize',14)
            ylabel('Eigenvector 2', 'FontName', 'Myriad Pro', 'FontSize',14)
            zlabel('Eigenvector 3', 'FontName', 'Myriad Pro', 'FontSize',14)
            xticklabels([])
            yticklabels([])
            zticklabels([])

            fg.Renderer = 'painters';
            fg.Position = [10 89 1304 1084];



            switch roiList{rr}
                case 'dACC'
                    set(gca, 'CameraPosition', [-0.3566   -1.2990    0.3254]);
                case 'PFCl'
                    set(gca, 'CameraPosition', [-0.6046   -1.1104    0.5491]);
                case 'SPL'
                    set(gca, 'CameraPosition', [-0.9856   -2.0389    0.7129]);
                case 'IPS'
                    set(gca, 'CameraPosition', [-0.7709   -2.0441    0.6028]);

            end

            %             ylim([-.1, .1])


        end


    end










end







%% MDS WITH RESP (3D IPS video)


addpath(genpath('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/pcm_toolbox'));
accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];
fefSel = [315, 251, 116, 58];
PFClSalSel = [103, 104, 105, 307, 308, 309];
PFClCtrlSel = [59, 60, 61, 62, 252, 253, 254, 268, 278, 279, 280];
accSalSel = [    85    86   106   286   287   295];
accCtrlSel = [    73    80    81   255   281   282   283];
motSel = [176, 375];
visaSel =  [172   173   174   175   176   177   178   179   180   181   182   369   370   371   372   373   374   375   376   377   378   379   380];
M1LeftSel = [149, 150, 151];
PrCdSel = [   116   130   131   315   329   330   331];
PFCdSel = [     5    20    21    22    23    32    33    45    46    58    72    92   115   203   217   218   219   220   229   238   251   267 306];


accLeftSel = [125,85,86,73,106,80,81];
accRightSel = [357,286,287,255,295,281,283,282];



load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/ScientificColourMaps7/roma/roma.mat')
cmap = roma;
cidx = round(linspace(1, 256, 8));


roiList = {'IPS'}
% roiList = {'dACC', 'PFCl', 'SPL', 'IPS'}
subnetworkList = {'ContA', 'ContB', 'ContC', 'DorsAttnA', 'DorsAttnB', 'SalVenAttnA', 'SalVenAttnB', 'SomMotA', 'SalVenAttnB', 'VisualA', 'VisualB', 'VisualC'}
networkList = {'Cont', 'DorsAttn', 'SalVenAttn', 'SomMot', 'Visual'}


if strcmp(analysisName, 'margResp_parcel')





        fg=figure;
 
        tiledlayout('flow', 'padding', 'compact', 'TileSpacing','compact')




    % higher dim
    for nn = 1

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




            % print stats
            M = mean(reshape(squeeze(mean(R(:,:,regSel,:),3)), [nCond^2, npt])');
            rM = reshape(M, [nCond, nCond]);




            [L,l] = pcm_classicalMDS(rM);

%                                     nexttile; imagesc(rM, [-.02, .02]); title(roiList{rr}); colormap(roma)
%                         nexttile; imagesc(L, [-.1, .1]); title(roiList{rr}); colormap(roma)



  
%             [L,l] = pcm_classicalMDS(rM,...
%                 'contrast', [...
%                 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0;...
%                 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1;...
%                 0.5, 0.1667, -0.1667, -0.5000, -0.5000, -0.1667, 0.1667, 0.5000, 0,0,0,0]');



            % combined



            mds = L;
            nexttile; hold on;


            % plot lines

            plot3(...
                mds(1:4,1),...
                mds(1:4,2),...
                mds(1:4,3), '-', 'color', cmap(cidx(1),:),...
                'LineWidth', 3 ...
                )

            plot3(...
                mds(5:8,1),...
                mds(5:8,2),...
                mds(5:8,3), '-', 'color', cmap(cidx(end),:),...
                'LineWidth', 3 ...
                )


            plot3(...
                [mean(mds([9,12],1)), mean(mds([10,11],1))],...
                [mean(mds([9,12],2)), mean(mds([10,11],2))],...
                [mean(mds([9,12],3)), mean(mds([10,11],3))],...
                '-k', 'LineWidth', 3 ...
                )


            for ii = 1:8
                plot3(mds(ii,1), mds(ii,2), mds(ii,3),  'o',...
                    'Color', cmap(cidx(ii),:), 'MarkerFaceColor', cmap(cidx(ii),:), 'MarkerSize', 18, 'MarkerEdgeColor','k');

                %                 plot3(mds(ii,1), zeros(size(mds(ii,1))) + .3, mds(ii,3),'o', 'MarkerSize', 10,'color',[0.8 0.8 0.8]);

            end

            for ii = 1:4
                plot3(mds(8+ii,1), mds(8+ii,2), mds(8+ii,3), 'd',...
                    'Color', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerFaceColor', cmap(cidx(2 + 2*(ii-1)),:), 'MarkerSize', 18, 'MarkerEdgeColor','k');

                %                 plot3(mds(8+ii,1), zeros(size(mds(8+ii,1))) + .3, mds(ii,3),'o', 'MarkerSize', 10,'color',[0.8 0.8 0.8]);

            end








            % shadow
            xl = xlim*1.1;
            yl = ylim*1.1;
            zl = zlim*1.1;

            shadowCol = [.75, .75, .75];

            zcol1 =  0*cmap(cidx(1),:) + .85*[.85, .85, .85];
            zcol2 =  0*cmap(cidx(end),:) + .85*[.85, .85, .85];
            zcol3 = 0*[0,0,0] + .85*[.85, .85, .85];


            col1 =  .15*cmap(cidx(1),:) + .85*[.9, .9, .9];
            col2 =  .15*cmap(cidx(end),:) + .85*[.9, .9, .9];
            col3 = .15*[0,0,0] + .85*[.9, .9, .9];


            plot3(...
                mds(1:4,1),...
                mds(1:4,2),...
                0*mds(1:4,3) + zl(1), 'o', ...
                'color', zcol1,...
                'MarkerSize', 18, ...
                'MarkerFaceColor', zcol1...
                )


            plot3(...
                mds(5:8,1),...
                mds(5:8,2),...
                0*mds(5:8,3) + zl(1), 'o', ...
                'color', zcol2,...
                'MarkerSize', 18, ...
                'MarkerFaceColor', zcol2...
                )


             plot3(...
                mds(9:10,1),...
                mds(9:10,2),...
                0*mds(9:10,3) + zl(1), 'd', ...
                'color', zcol1,...
                'MarkerSize', 18, ...
                'MarkerFaceColor', zcol1...
                )

             plot3(...
                 mds(11:12,1),...
                 mds(11:12,2),...
                 0*mds(11:12,3) + zl(1), 'd', ...
                 'color', zcol2,...
                 'MarkerSize', 18, ...
                 'MarkerFaceColor', zcol2...
                 )


            % z project
            plot3(...
                mds(1:4,1),...
                mds(1:4,2),...
                0*mds(1:4,3) + zl(1), '-', ...
                'LineWidth', 3,...
                'color', zcol1,...
                'MarkerFaceColor', zcol1...
                )

            plot3(...
                mds(5:8,1),...
                mds(5:8,2),...
                0*mds(5:8,3) + zl(1), '-', ...
                'LineWidth', 3,...
                'color', zcol2,...
                'MarkerFaceColor', zcol2...
                )


            plot3(...
                [mean(mds([9,12],1)), mean(mds([10,11],1))],...
                [mean(mds([9,12],2)), mean(mds([10,11],2))],...
                [0*mean(mds([9,12],3)) + zl(1), 0*mean(mds([10,11],3)) + zl(1)], '-', ...
                'LineWidth', 3,...
                'color', zcol3,...
                'MarkerFaceColor', zcol3...
                )



            % y project
%             plot3(...
%                 mds(1:4,1),...
%                 0*mds(1:4,2) + yl(2),...
%                 mds(1:4,3), '-', ...
%                 'LineWidth', 2,...
%                 'color', col1,...
%                 'MarkerFaceColor', col1...
%                 )
% 
%             plot3(...
%                 mds(5:8,1),...
%                 0*mds(5:8,2) + yl(2),...
%                 mds(5:8,3), '-', ...
%                 'LineWidth', 2,...
%                 'color', col2,...
%                 'MarkerFaceColor', col2...
%                 )
% 
% 
%             plot3(...
%                 [mean(mds([9,12],1)), mean(mds([10,11],1))],...
%                 [0*mean(mds([9,12],2)) + yl(2), 0*mean(mds([10,11],2))+ yl(2)],...
%                 [mean(mds([9,12],3)), mean(mds([10,11],3))], '-', ...
%                 'LineWidth', 2,...
%                 'color', col3,...
%                 'MarkerFaceColor', col3...
%                 )
% 
% 
% 
% 
% 
%             % x project
%             plot3(...
%                 0*mds(1:4,1) + xl(2),...
%                 mds(1:4,2),...
%                 mds(1:4,3), '-', ...
%                 'LineWidth', 2,...
%                 'color', col1,...
%                 'MarkerFaceColor', col1...
%                 )
% 
%             plot3(...
%                 0*mds(5:8,1) + xl(2),...
%                 mds(5:8,2),...
%                 mds(5:8,3), '-', ...
%                 'LineWidth', 2,...
%                 'color', col2,...
%                 'MarkerFaceColor', col2...
%                 )
% 
% 
%             plot3(...
%                 [0*mean(mds([9,12],1)) + xl(2), 0*mean(mds([10,11],1))+ xl(2)],...
%                 [mean(mds([9,12],2)), mean(mds([10,11],2))],...
%                 [mean(mds([9,12],3)), mean(mds([10,11],3))], '-', ...
%                 'LineWidth', 2,...
%                 'color', col3,...
%                 'MarkerFaceColor', col3...
%                 )




            %             title([regionList{rr}])
            set(gca, 'TickDir', 'None', 'LineWidth', 1.5)
            axis('square')
            grid on;

            xlim(xl)
            ylim(yl)
            zlim(zl)
            xlabel('Eigenvector 1', 'FontName', 'Myriad Pro', 'FontSize',14)
            ylabel('Eigenvector 2', 'FontName', 'Myriad Pro', 'FontSize',14)
            zlabel('Eigenvector 3', 'FontName', 'Myriad Pro', 'FontSize',14)
            xticklabels([])
            yticklabels([])
            zticklabels([])

            fg.Renderer = 'painters';
            fg.Position = [10 89 1304 1084];


            text1=  text(mds(2,1), mds(2,2)+.075, mds(2,3), 'Target Coherence (Right)', 'color', cmap(cidx(1),:), 'FontSize', 18,'FontName', 'Myriad Pro');
            text2= text(mds(6,1), mds(6,2)+.075, mds(6,3), 'Target Coherence (Left)', 'color', cmap(cidx(end),:), 'FontSize', 18,'FontName', 'Myriad Pro');
            text3= text((mds(9,1) + mds(10,1))/2, (mds(9,2) + mds(10,2))/2 +.1, (mds(9,3) + mds(10,3))/2, 'Distractor Coherence', 'color', 'k', 'FontSize', 18,'FontName', 'Myriad Pro');



            init_pos = [-0.7709   -1.5    0.6028];

            set(gca, 'CameraPosition', init_pos);
%             exportgraphics(gcf,'testAnimated.gif');

            nframe = 30;
            fxx = linspace(init_pos(1), 0, nframe);
            fyy = linspace(init_pos(2), -5, nframe);
            fzz = linspace(init_pos(3), 10, nframe);


            for ii = 1:15

                set(gca, 'CameraPosition', [fxx(1)   fyy(end)    fzz(1)]);
                exportgraphics(gcf,'testAnimated.gif','Append',true);
            end

            delete(text1)
            delete(text2)
            delete(text3)


            for jj = 1:nframe
               
                set(gca, 'CameraPosition', [fxx(jj)   fyy(end)    fzz(jj)]);
                exportgraphics(gcf,'testAnimated.gif','Append',true);
            end

            for jj = nframe:-1:1
                set(gca, 'CameraPosition', [fxx(jj)   fyy(end)    fzz(jj)]);
                exportgraphics(gcf,'testAnimated.gif','Append',true);
            end

%             for ii = nframe:-1:1
% 
%                 set(gca, 'CameraPosition', [fxx(1)   fyy(ii)    fzz(1)]);
%                 exportgraphics(gcf,'testAnimated.gif','Append',true);
%             end




            %             ylim([-.1, .1])


        end


    end










end












%% compare between tasks


if strcmp(analysisName, 'featureBlkTask_parcel') | strcmp(analysisName, 'perfBlkTask_parcel') | strcmp(analysisName, 'perfFullBlkTask_parcel')



    mx_sims = 10000;
    hb_alpha = .05;

    alpha = .001;
    condSel = @(x) find(strcmp(condLabel, x));

    mkdir(fullfile(save_dir, 'taskComp'));
    delete(fullfile(save_dir, 'taskComp', '*.nii'))


    for mm = 3 % for each metric

        % G reliability
        for vv = 1:length(corrName)

            % get stats
            switch mm
                case 1


                case 2


                case 3

                    dists = squeeze(R(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{1}),:,:))' - squeeze(R(condSel(corrPairs{vv}{2}),condSel(corrPairs{vv}{2}),:,:))';

                    [~,pval,~,stats] = ttest(squeeze(R(condSel(corrPairs{vv}{1}),condSel(corrPairs{vv}{1}),:,:))', squeeze(R(condSel(corrPairs{vv}{2}),condSel(corrPairs{vv}{2}),:,:))');
                    fnameT = fullfile(save_dir,'taskComp', sprintf('T_%s.nii', corrName{vv}));
                    fnameP = fullfile(save_dir,'taskComp', sprintf('phb_%s.nii', corrName{vv}));
                    fnameBin = fullfile(save_dir,'taskComp', sprintf('mx_bin_%s.nii', corrName{vv}));

            end

            %             [psort, idxsort] = sort(pval);
            %             hb_pval = zeros(1, length(pval));
            %             hb_pval(idxsort) = psort .*[length(psort):-1:1];
            %             %             hb_pval = pval*17;
            %             hb_pval = min(hb_pval, 1);


            % get correction

            [psort, idxsort] = sort(pval);
            hb_pval = zeros(1, length(pval));
            corP = psort .*[length(psort):-1:1];
            corP(find(corP>.05, 1):end) = 1;
            hb_pval(idxsort) = corP;





            % randomization test
            perm_mx = reshape(datasample([-1,1], numel(dists)*mx_sims), [size(dists,1), size(dists,2), mx_sims]);
            null_R = repmat(dists, [1,1,mx_sims]) .* perm_mx;

            max_null_R = max(squeeze(mean(null_R) ./ std(null_R)));
            min_null_R = min(squeeze(mean(null_R) ./ std(null_R)));
            true_R = mean(dists)./std(dists);
            if isempty(max_null_R)
                max_null_R = 0;
                min_null_R = 0;
            end





            %             [~, ~, ~, hb_pval]=fdr_bh(pval);

            % print to parcellation
            parcelBrain = maskImg;
            for rr = 1:nParcel
                parcelBrain(parcelBrain == rr) = stats.tstat(rr);
            end

            % write
            VoOut      = struct(...
                'fname',    fnameT,...
                'dim',      maskVo.dim,...
                'dt',       [spm_type('float32') spm_platform('bigend')],...
                'mat',      maskVo.mat .* [-1 1 1 1]',...
                'n',        [1 1],...
                'descrip',  'task contrast');

            spm_write_vol(VoOut, parcelBrain);




            % print to parcellation
            parcelBrain = maskImg;
            for rr = 1:nParcel
                parcelBrain(parcelBrain == rr) = -log10(hb_pval(rr));
            end

            % write
            VoOut      = struct(...
                'fname',    fnameP,...
                'dim',      maskVo.dim,...
                'dt',       [spm_type('float32') spm_platform('bigend')],...
                'mat',      maskVo.mat .* [-1 1 1 1]',...
                'n',        [1 1],...
                'descrip',  'pvalue');

            spm_write_vol(VoOut, parcelBrain);



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
                'fname',    fnameBin,...
                'dim',      maskVo.dim,...
                'dt',       [spm_type('float32') spm_platform('bigend')],...
                'mat',      maskVo.mat .* [-1 1 1 1]',...
                'n',        [1 1],...
                'descrip',  'pvalue');

            spm_write_vol(VoOut, parcelBrain);





        end
    end




end









%% compare between tasks (REGION)




% =========


if strcmp(analysisName, 'featureBlkTask_parcel') | strcmp(analysisName, 'perfBlkTask_parcel') | strcmp(analysisName, 'perfFullBlkTask_parcel')



    accSel = [357,286,287,255,295,281,283,282,125,85,86,73,106,80,81];
    fefSel = [315, 251, 116, 58];
    PFClSalSel = [103, 104, 105, 307, 308, 309]
    PFClCtrlSel = [59, 60, 61, 62, 252, 253, 254, 268, 278, 279, 280]
    accSalSel = [    85    86   106   286   287   295]
    accCtrlSel = [    73    80    81   255   281   282   283]


    % roiList = {'ExStrSup', 'SPL', 'IPS', 'IPL', 'PFCl', 'dACC'}
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
                case 'FEF'
                    regSel = ismember(1:400, fefSel);
                case 'PFClSal'
                    regSel = ismember(1:400, PFClSalSel);
                case 'PFClCtrl'
                    regSel = ismember(1:400, PFClCtrlSel);
                case 'accSal'
                    regSel = ismember(1:400, accSalSel);
                case 'accCtrl'
                    regSel = ismember(1:400, accCtrlSel);
                otherwise
                    regSel = cellfun(@any, regexpi(regionTbl, regionList{rr}));
            end


            col_ms = [...
                squeeze(mean(R(condSel('colRespTarg'),condSel('colRespTarg'),regSel,:),3)),...
                squeeze(mean(R(condSel('colRespDist'),condSel('colRespDist'),regSel,:),3)),...
                squeeze(mean(R(condSel('colCohTarg'),condSel('colCohTarg'),regSel,:),3)),...
                squeeze(mean(R(condSel('colCohDist'),condSel('colCohDist'),regSel,:),3)),...
                squeeze(mean(R(condSel('colCong'),condSel('colCong'),regSel,:),3)),...
                ];

            mot_ms = [...
                squeeze(mean(R(condSel('motRespTarg'),condSel('motRespTarg'),regSel,:),3)),...
                squeeze(mean(R(condSel('motRespDist'),condSel('motRespDist'),regSel,:),3)),...
                squeeze(mean(R(condSel('motCohTarg'),condSel('motCohTarg'),regSel,:),3)),...
                squeeze(mean(R(condSel('motCohDist'),condSel('motCohDist'),regSel,:),3)),...
                squeeze(mean(R(condSel('motCong'),condSel('motCong'),regSel,:),3)),...
                ];


            disp(regionList{rr})
            disp('=============')

            [~,col_pval,~,col_stats] = ttest(col_ms);
            [~,mot_pval,~,mot_stats] = ttest(mot_ms);
            [~,comp_pval,~,comp_stats] = ttest(col_ms, mot_ms);
            [~,resp_pval,~,resp_stats] = ttest((col_ms(:,1) + mot_ms(:,1)), (col_ms(:,2) + mot_ms(:,2)), 'Tail', 'right')
            [~,col_resp_pval,col_resp_CI,col_resp_stats] = ttest(col_ms(:,1), (col_ms(:,2)), 'Tail', 'right')
            col_resp_D = col_resp_stats.tstat./sqrt(col_resp_stats.df+1)
            [~,mot_resp_pval,mot_resp_CI,mot_resp_stats] = ttest((mot_ms(:,1)), (mot_ms(:,2)), 'Tail', 'right')
            mot_resp_D = mot_resp_stats.tstat./sqrt(mot_resp_stats.df+1)

            disp('=============')

            nexttile; hold on;
            errorbar([1:5]-.15, mean(col_ms), std(col_ms)./sqrt(npt),...
                'o', 'Color', taskCols{1}, 'MarkerFaceColor', taskCols{1}, 'MarkerSize', 10, 'LineWidth', 2);
            errorbar([1:5]+.1, mean(mot_ms), std(mot_ms)./sqrt(npt),...
                'o', 'Color', taskCols{2}, 'MarkerFaceColor', taskCols{2},'MarkerSize', 10, 'LineWidth', 2);

            title(regionList{rr})
            xlim([.5, 5.5])
            xticks(1:5)
            xticklabels(condLabel(1:end-1))
            yline(0)
            ylim([-.02, .10])
            set(gca, 'TickDir', 'out', 'LineWidth', 1)





        end






    end

end

















%% compare SPL-IPS FC


compare_name = 'TargCohComp'
feature_idx = 3;




analysisName = 'FC_IPS_parcel' %

% get paths
save_dir = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s',analysisFolder, analysisName)
data_dir = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s/fit-results',analysisFolder, analysisName)

dr = dir(fullfile(data_dir, '*.mat'));



% load data
[G,C,R1,Rw] = deal([]);

for ii = 1:length(dr)

    f = load(fullfile(dr(ii).folder, dr(ii).name));

    G = cat(4, G, f.G);
    C = cat(4, C, f.C);
    R1 = cat(4, R1, f.R);


end



analysisName = 'FC_SPL_parcel' %

% get paths
data_dir = sprintf('/Users/hr0283/Dropbox (Brown)/RDM_fmri_results/%s/%s/fit-results',analysisFolder, analysisName)

dr = dir(fullfile(data_dir, '*.mat'));



[G,C,R2,Rw] = deal([]);

for ii = 1:length(dr)

    f = load(fullfile(dr(ii).folder, dr(ii).name));

    G = cat(4, G, f.G);
    C = cat(4, C, f.C);
    R2 = cat(4, R2, f.R);


end





mx_sims = 10000;
hb_alpha = .05;

alpha = .001;
condSel = @(x) find(strcmp(condLabel, x));


mkdir(fullfile(save_dir, compare_name));
delete(fullfile(save_dir, compare_name, '*.nii'))


vv=1;


dists = squeeze(R2(feature_idx,end,:,:))' - squeeze(R1(feature_idx,end,:,:))';

[~,pval,~,stats] = ttest(squeeze(R2(feature_idx,end,:,:))', squeeze(R1(feature_idx,end,:,:))');
fnameT = fullfile(save_dir,compare_name, sprintf('T_%s.nii',compare_name));
fnameP = fullfile(save_dir,compare_name, sprintf('phb_%s.nii', compare_name));
fnameBin = fullfile(save_dir,compare_name, sprintf('mx_bin_%s.nii', compare_name));



% get correction

[psort, idxsort] = sort(pval);
hb_pval = zeros(1, length(pval));
corP = psort .*[length(psort):-1:1];
corP(find(corP>.05, 1):end) = 1;
hb_pval(idxsort) = corP;





% randomization test
perm_mx = reshape(datasample([-1,1], numel(dists)*mx_sims), [size(dists,1), size(dists,2), mx_sims]);
null_R = repmat(dists, [1,1,mx_sims]) .* perm_mx;

max_null_R = max(squeeze(mean(null_R) ./ std(null_R)));
min_null_R = min(squeeze(mean(null_R) ./ std(null_R)));
true_R = mean(dists)./std(dists);
if isempty(max_null_R)
    max_null_R = 0;
    min_null_R = 0;
end





%             [~, ~, ~, hb_pval]=fdr_bh(pval);

% print to parcellation
parcelBrain = maskImg;
for rr = 1:nParcel
    parcelBrain(parcelBrain == rr) = stats.tstat(rr);
end

% write
VoOut      = struct(...
    'fname',    fnameT,...
    'dim',      maskVo.dim,...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      maskVo.mat .* [-1 1 1 1]',...
    'n',        [1 1],...
    'descrip',  'task contrast');

spm_write_vol(VoOut, parcelBrain);




% print to parcellation
parcelBrain = maskImg;
for rr = 1:nParcel
    parcelBrain(parcelBrain == rr) = -log10(hb_pval(rr));
end

% write
VoOut      = struct(...
    'fname',    fnameP,...
    'dim',      maskVo.dim,...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      maskVo.mat .* [-1 1 1 1]',...
    'n',        [1 1],...
    'descrip',  'pvalue');

spm_write_vol(VoOut, parcelBrain);



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
    'fname',    fnameBin,...
    'dim',      maskVo.dim,...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      maskVo.mat .* [-1 1 1 1]',...
    'n',        [1 1],...
    'descrip',  'pvalue');

spm_write_vol(VoOut, parcelBrain);










































