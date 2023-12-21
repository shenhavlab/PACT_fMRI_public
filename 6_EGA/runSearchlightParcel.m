function runSearchlightParcel(parcel,varargin)
% runSearchlightLDC(searchLight,varargin)
% Wrapper for the main searchlight function rsa_runSearchlight
% This version calculates the LDC from an SPM first-level analysis
% New Version takes into account the first-level design matrix
%
% INPUT:
%     searchLight:  a) Name of the precalculated searchlight file
%                   b) Searchlight structure itself (see rsa_defineSearchlight.m)
%
% VARARGIN: either list of 'optionname',value,'optioname2',value2..
%           or structure with Opt.optioname = value
%   rootPath:       Path to be pre-pended to all path
%   spmDir:         Path (relative to root) where the SPM.mat resides and
%                   were results are written
%   conditionVec:   Vector of condition labels for the task-related
%                   regressors in the SPM matrix, if empty, the routine
%                   assumes that all betas are condition1.. conditionK
%   conditionLabels:Cell array of condition labels. If given, the routine
%                   attempts to extract the conditions and partitions from
%                   the SPM structure
%   imageDataDir:   Alternative directory to find the raw imaging data, if
%                   it different from the one that is specified in the SPM
%                   structure. This is especially useful if the SPM was
%                   estimated in one directory and then later moved.
% (C) Joern Diedrichsen 2015

import rsa.spm.*
import rsa.util.*

disp(' =============== runSearchlightParcel =============== ')


% User optional parameters
Opt.rootPath        = [];
Opt.spmDir          = [];     %  Directory that contains the first-level analysis
Opt.analysisName    = [];     %  Name of the analysis
Opt.conditionVec    = [];     %  Vector, indicating which of the betas belong to which condition
%  Use 0 for betas that will not be included
Opt.conditionLabels = {};     %  Condition labels - used if conditionVec is not given
Opt.partition       = [];     %  Indicator of partition
Opt.imageDataDir    = [];   %  Directory where the pre-processed raw data resides (if different from what is specified in the SPM-structure)
Opt.saveSigma       = 'none'; %  Determine whether to save the variance / covariance of the beta estimates as well as the distances

Opt = rsa.getUserOptions(varargin,Opt);

workingDir = fullfile(Opt.rootPath, Opt.spmDir);

% load the parcellation
parcVol  = flipud(spm_read_vols(parcel.Vo));

% Load the SPM structure
s=load(fullfile(workingDir,'SPM.mat'));
SPM = s.SPM

% report collinearity
collintest(SPM.xX.xKXs.X(:, find(Opt.conditionVec~=0)))


% Determine the condition and the partition labels
nBetas = SPM.xX.iB - Opt.nBlk; % Use only task related regressors (no intercepts)
if (isempty(Opt.conditionVec))
    if (isempty(Opt.conditionLabels))
        error('Either Opt.conditionVec or Opt.conditionLabels needs to be set');
    end;
    [Opt.conditionVec,Opt.partition]=rsa.getSPMconditionVec(SPM,Opt.conditionLabels);
else
    if (isempty(Opt.partition))
        Opt.partition=zeros(nBetas,1);
        for i=1:length(SPM.Sess)
            Opt.partition(SPM.Sess(i).col)=i;
        end;
    end;
end;

Opt.partition;

nConditions = max(Opt.conditionVec);

% Now define the input files
if (isempty(Opt.imageDataDir))
    inFiles = SPM.xY.VY;
else
    numFiles = length(SPM.xY.VY);
    for i=1:numFiles
        [dir,filename,extension,number]=spm_fileparts(SPM.xY.P(i,:));
        inFileNames{i} = fullfile(Opt.rootPath,Opt.imageDataDir,sprintf('%s%s%s',filename,extension,number));
    end;
    inFiles = spm_vol(char(inFileNames));
end;

if (isstruct(inFiles))
    VolIn=inFiles;
else
    VolIn=spm_vol(inFiles);
end;



% organize stuff & get similarities
partition       = Opt.partition;
conditionVec    = Opt.conditionVec;
voDim           = parcel.Vo.dim;
nParcel         = height(parcel.info);

[G,C,T,R,S,absR] = deal(nan(nConditions, nConditions, nParcel));

fprintf('\nstarting fit ---------\n')

tic;
for rr = 1:nParcel

    % get idx
    [I,J,K] = ind2sub(voDim, find(parcVol == rr));

    % get data
    nVx(rr) = length(I);
    Y = nan(nVx(rr), length(VolIn));
    for i=1:length(VolIn)
        Y(:,i)=spm_sample_vol(VolIn(i),double(I),double(J),double(K),0);
    end

    % get similarities
    [G(:,:,rr), C(:,:,rr), T(:,:,rr), R(:,:,rr), S(:,:,rr), absR(:,:,rr)]  = calculateCov(Y', SPM, partition, conditionVec);

    % display progress
    if mod(rr, round(nParcel/5)) == 0
        fprintf('--- parcel %d/%d. dur: %.2gm\n', rr, nParcel, toc/60)
        tic;
    end

end
fprintf('\finished fit ---------\n')




% get info
[subNetwork, network, region] = deal({});
for rr = 1:height(parcel.info)
    name = strsplit(parcel.info.Var5{rr}, '_');
    subNetwork{rr,1} = name{3};
    network{rr,1}    = name{3}(1:end-1);
    region{rr,1}     = name{4};
end

info = cell2table([subNetwork, network, region], 'VariableNames', {'subNetwork', 'network', 'region'});


% === SAVE
save(Opt.saveFn, 'G', 'C', 'T', 'R', 'absR', 'S', 'nVx', 'info', 'parcel', 'Opt')



end





%% plug-in function for rsa_runSearchlight
function [G, C, T, R, S, absR] = calculateCov(Y, SPM, partition, conditionVec)

% get (whitened) similarity
U_hat   = rsa.spm.noiseNormalizeBeta(Y,SPM); % Get noise normalised betas
[G, C, T, R, S, absR] = pcm_estGCrossval(U_hat, partition, conditionVec, 'X', SPM.xX.xKXs.X); % get crossval G


end








