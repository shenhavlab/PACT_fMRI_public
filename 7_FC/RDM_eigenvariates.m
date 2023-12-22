function extract_eigenvariates(root_dir, spm_dir, ptNum, varargin)
%% ROI analysis
% Harrison Ritz 2022




%% === set variables
addpath(genpath(spm_dir)); % add spm12 folder with bug-squashed rwls
spm('defaults', 'fmri')

% set default RAM & use RAM for analysis
spm_get_defaults('maxmem', 128 * 2^30)
spm_get_defaults('resmem',true)
spm_get_defaults('cmdline',true)



%% load variables
if length(varargin) >= 1 && ~isempty(varargin{1})
    name = varargin{1};
else
    name        = 'feature';
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

pt_dir      = sprintf('%s/spm-data/sub-%d', root_dir, ptNum)
spm_dir     = sprintf('%s/level 1/%s', pt_dir, analysis)
mask_dir    = sprintf('%s/RDM_fmri_scripts/masks', root_dir)
save_dir    = sprintf('%s/spm-data/FC/data/', root_dir)


% make directories
mkdir(save_dir); % save results




%% regions

regList = {...
    [24, 47, 59, 60, 61, 62, 103, 104, 105, 239, 240, 252, 253, 254, 268, 278, 279, 280, 307, 308, 309],... % PFCl
    [53, 54, 55, 67, 111, 112, 113, 114, 246, 247, 248, 249, 262, 313, 314],... IPS
    }


regName = {...
    [name, '_PFCl',]... % PFCl
    [name, '_IPS',]... % HG
    }



psych = {'targCoh', 'distCoh', 'cong'}


%% get data
% SPM mat
load(fullfile(spm_dir, 'SPM.mat'))

% parcellation
parcel.Vo   = spm_vol(fullfile(mask_dir, 'Schaefer2018_400Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii'));
parcel.info = readtable(fullfile(mask_dir, 'Schaefer2018_400Parcels_Kong2022_17Networks_order.lut.txt'));
parcel
parcVol  = flipud(spm_read_vols(parcel.Vo)); % fuck, need to flip parcel

% image loc
inFiles = SPM.xY.VY;
if (isstruct(inFiles))
    VolIn=inFiles;
else
    VolIn=spm_vol(inFiles);
end



%% get timeseries [ ======= built on SPM PPI code  ======= ]

nBlk = (6 - ismember(ptNum, [8010, 8019]));
runMx = [];
for bb = 1:nBlk
    runMx = blkdiag(runMx, ones(392,1));
end
runMx = logical(runMx);

pc1 = double(runMx(:,1))*0;
pc2 = double(runMx(:,1))*0;

for rr = 1:length(regList)

    % ======= get beta & multivariate noise norm ============================

    [I,J,K] = ind2sub(parcel.Vo.dim, find(ismember(parcVol, regList{rr})));

    % get data
    Y = nan(length(I), length(VolIn));

    for vv = 1:length(VolIn)
        Y(:,vv)=spm_sample_vol(VolIn(vv), double(I), double(J), double(K), 0);
    end


    % get res
    KWY     = spm_filter(SPM.xX.K, SPM.xX.W*Y');        %%% filter out low-frequence trends in Y
    res     = spm_sp('r', SPM.xX.xKXs, KWY);            %%% residuals: res  = Y - X*beta

    szRes = size(res)
    res = res(:, all(isfinite(res),1)); % remove bad voxels
    szRes = size(res)

    % get eigenvariate
    for bb = 1:nBlk

        y = res(runMx(:,bb),:);

        [m,n]   = size(y);

        % PC1
        if m > n
            [v,s,v] = svd(y'*y);
            s       = diag(s);
            v       = v(:,1);
            u       = y*v/sqrt(s(1));
        else
            [u,s,u] = svd(y*y');
            s       = diag(s);
            u       = u(:,1);
            v       = y'*u/sqrt(s(1));
        end
        d       = sign(sum(v));
        u       = u*d;
        v       = v*d;
        pc1(runMx(:,bb)) = u*sqrt(s(1)/n);

        % PC2
        if m > n
            [v,s,v] = svd(y'*y);
            s       = diag(s);
            v       = v(:,2);
            u       = y*v/sqrt(s(2));
        else
            [u,s,u] = svd(y*y');
            s       = diag(s);
            u       = u(:,2);
            v       = y'*u/sqrt(s(2));
        end
        d       = sign(sum(v));
        u       = u*d;
        v       = v*d;
        pc2(runMx(:,bb)) = u*sqrt(s(2)/n);




    end

    szY = size(pc1)



    % do PPI ==================================================================

    % Setup variables
    %--------------------------------------------------------------------------
    Sess = SPM.Sess;
    RT      = SPM.xY.RT;
    dt      = SPM.xBF.dt;
    NT      = round(RT/dt);
    fMRI_T0 = SPM.xBF.T0;
    N       = length(pc1);
    k       = 1:NT:N*NT;                       % microtime to scan time indices

    % Create basis functions and hrf in scan time and microtime
    %--------------------------------------------------------------------------
    hrf = spm_hrf(dt);

    % Create convolved explanatory {Hxb} variables in scan time
    %--------------------------------------------------------------------------
    xb  = spm_dctmtx(N*NT + 128,N);
    Hxb = zeros(N,N);
    for i = 1:N
        Hx       = conv(xb(:,i),hrf);
        Hxb(:,i) = Hx(k + 128);
    end
    xb = xb(129:end,:);

    % Specify covariance components; assume neuronal response is white
    % treating confounds as fixed effects
    %--------------------------------------------------------------------------
    Q = speye(N,N)*N/trace(Hxb'*Hxb);
    % Q = blkdiag(Q, speye(M,M)*1e6  );

    % Get whitening matrix (NB: confounds have already been whitened)
    %--------------------------------------------------------------------------
    W = SPM.xX.W(Sess.row,Sess.row);

    % Create structure for spm_PEB
    %--------------------------------------------------------------------------
    clear P
    P{1}.X = W*Hxb;         % Design matrix for lowest level
    P{1}.C = speye(N,N)/4;  % i.i.d assumptions
    P{2}.X = sparse(N,1);   % Design matrix for parameters (0's)
    P{2}.C = Q;

    % do deconv with spm_PEB
    C  = spm_PEB(pc1,P);
    xn = xb*C{2}.E(1:N);
    xn = spm_detrend(xn);

    % calc PPI
    clear psy
    for pp = 1:length(psych)

        % Setup psychological variable from inputs and contrast weights
        %----------------------------------------------------------------------
        PSY = full(Sess.U(1).u(33:end, ~cellfun(@isempty, regexpi(Sess.U(1).name, psych{pp}))));

        % Multiply psychological variable by neural signal
        %----------------------------------------------------------------------
        PSYxn = PSY.*xn;

        % Convolve, convert to scan time, and account for slice timing shift
        %----------------------------------------------------------------------
        for ii = 1:size(PSYxn,2)
            blk_ppi = conv(PSYxn(:,ii),hrf);
            blk_ppi = blk_ppi((k-1) + fMRI_T0);
            psy(pp).(psych{pp})(:,ii) = blk_ppi;
        end

    end

    % save
    savename = fullfile(save_dir, sprintf('%d_%s.mat', ptNum, regName{rr}))
    save(savename, 'pc1', 'pc2', 'psy')



end



