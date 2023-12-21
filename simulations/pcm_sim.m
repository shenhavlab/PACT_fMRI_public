%%
clear;clc

addpath(genpath('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/pcm_toolbox'))

%% parameters
nTrials = 150;
D.numPart = 6;
D.numVox  = 250;
numSim = 100;
signal = 1;
noise  = 1;

%% model

M.type       = 'feature';
M.numGparams = 2;
M.Ac(:,:,1) = diag([1, 0]);
M.Ac(:,:,2) = diag([0, 1]);
M.theta0=ones(2,1);                        % Starting values: could be closer, but converges anyways 


%% design matrix
collin = .3;
X = mvnrnd([0,0], [1, collin; collin, 1], nTrials);

%% generate
[Y, partVec, condVec] = pcm_generateData(M, M.theta0, D, numSim, signal, noise, 'design', X);

%% fit PCM
tic
[Z,B,X,YY,S,N,P,G_hat,noise0,run0]=...
    pcm_setUpFit(Y,partVec,condVec,'runEffect','random');
toc

squeeze(mean(G_hat,3))


[BF,T]=bf_ttest(squeeze(G_hat(1,2,:)))
logBF = log10(BF)
























%% within sample ================= PAPER SIMULATION

bam = load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/ScientificColourMaps7/bam/bam.mat')

%% parameters
nTrials = 150;
D.numPart = 1;
D.numVox  = 250;
numSim = 10000;
signal = 1;
noise  = 1;

%% model

enc2 = 0;
rel = .5;

clear M
M.type       = 'feature';
M.numGparams = 2;
M.Ac(:,:,1) = [1 0 rel 0; 0 0 0 0; rel 0 1 0; 0 0 0 0];
M.Ac(:,:,2) = [0 0 0 0; 0 enc2 0 enc2*rel;  0 0 0 0; 0 enc2*rel 0 enc2];
M.theta0=ones(2,1);                        % Starting values: could be closer, but converges anyways 


%% design matrix
collin = .5;
X = blkdiag(mvnrnd([0,0], [1, collin; collin, 1], nTrials), mvnrnd([0,0], [1, collin; collin, 1], nTrials));

%% generate
[Y, partVec, condVec] = pcm_generateData(M, M.theta0, D, numSim, signal, noise, 'design', X);


%% get similarity
R = zeros(4,4,numSim);
for ii = 1:length(Y)

   B = condVec \ Y{ii};
   R(:,:,ii) = corr(B');

end


%% plot

figure;
colormap(bam.bam)

nexttile;
imagesc(condVec, [-3,3])
colorbar;
axis('square')

nexttile;
imagesc(corr(condVec), [-1,1])
colorbar;
axis('square')

nexttile;
imagesc(M.Ac(:,:,1), [-1,1])
colorbar;
axis('square')

nexttile;
imagesc(M.Ac(:,:,2), [-1,1])
colorbar;
axis('square')

nexttile;
imagesc((condVec \ Y{1})', [-3,3])
colorbar;
axis('square')

nexttile;
imagesc(mean(R,3)./std(R,[],3), [-1 1])
axis('square')























%% absR simulation, combined =================================

bam = load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/ScientificColourMaps7/bam/bam.mat')

%% parameters
nTrials = 1000;
D.numPart = 6;
D.numVox  = 250;
numSim = 1000;
signal = 1;
noise  = 1;

%% model

M.type       = 'feature';
M.numGparams = 2;
M.Ac(:,:,1) = diag([1, 0]);
M.Ac(:,:,2) = diag([0, 1]);
M.theta0=ones(2,1);                        % Starting values: could be closer, but converges anyways 


%% design matrix
collin = 0;
X = mvnrnd([0,0], [1, collin; collin, 1], nTrials);

%% generate
[Y, partVec, condVec] = pcm_generateData(M, M.theta0, D, numSim, signal, noise, 'design', X);


%% get similarity
R = zeros(2,2,numSim);
for ii = 1:length(Y)

   B = condVec \ Y{ii};
   R(:,:,ii) = corr(B');
   R(:,:,ii) = corr(abs(B)');

end


%% plot

figure;
colormap(bam.bam)

nexttile;
imagesc(condVec, [-3,3])
colorbar;
axis('square')

nexttile;
imagesc(corr(condVec), [-1,1])
colorbar;
axis('square')

nexttile;
imagesc(M.Ac(:,:,1), [-1,1])
colorbar;
axis('square')

nexttile;
imagesc(M.Ac(:,:,2), [-1,1])
colorbar;
axis('square')

nexttile;
imagesc((condVec \ Y{1})', [-3,3])
colorbar;
axis('square')

nexttile;
imagesc(mean(R,3)./std(R,[],3), [-1 1])
axis('square')










%% absR simulation, combined =================================

bam = load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/ScientificColourMaps7/bam/bam.mat')
roma = load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/ScientificColourMaps7/roma/roma.mat')
broc = load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/ScientificColourMaps7/broc/broc.mat')

%% parameters
nTrials = 150;

D_sep.numPart = 6;
D_sep.numVox  = 125;

D_comb.numPart = 6;
D_comb.numVox  = 250;


numSim = 10000;
signal = 1;
noise  = 1;

%% design matrix
collin = 0;
X = mvnrnd([0,0], [1, collin; collin, 1], nTrials);

%% model

M.type       = 'feature';
M.numGparams = 1;
M.Ac = diag([1, 0]);
M.theta0=ones(1,1);                        % Starting values: could be closer, but converges anyways 

[Y1, partVec, condVec] = pcm_generateData(M, M.theta0, D_sep, numSim, signal, noise, 'design', X);


M.type       = 'feature';
M.numGparams = 1;
M.Ac = diag([0, 1]);
M.theta0=ones(1,1);                        % Starting values: could be closer, but converges anyways 

[Y2, partVec, condVec] = pcm_generateData(M, M.theta0, D_sep, numSim, signal, noise, 'design', X);



M.type       = 'feature';
M.numGparams = 1;
M.Ac(:,:,1) = diag([1, 0]);
M.Ac(:,:,2) = diag([0, 1]);
M.theta0=ones(1,1);                        % Starting values: could be closer, but converges anyways 

[Ycomb, partVec, condVec] = pcm_generateData(M, M.theta0, D_comb, numSim, signal, noise, 'design', X);


%% get similarity
[R_sep, absR_sep,...
    R_isol_sep, absR_isol_sep,...
    R_comb, absR_comb,...
    R_isol, absR_isol] = deal(zeros(2,2,numSim));


tic
for ii = 1:length(Y1)

   B = condVec \ [Y1{ii}, Y2{ii}];
   R_sep(:,:,ii) = corr(B');
   absR_sep(:,:,ii) = corr(abs(B)');


   y10 = Y1{ii}; y10(:, 1:50) = 0;
   y20 = Y2{ii}; y20(:, 1:50) = 0;
   B = condVec \ [y10, y20];
   R_isol_sep(:,:,ii) = corr(B');
   absR_isol_sep(:,:,ii) = corr(abs(B)');


   B = condVec \ [Ycomb{ii}];
   R_comb(:,:,ii) = corr(B');
   absR_comb(:,:,ii) = corr(abs(B)');


   B = condVec \ [Ycomb{ii}];
   B(:,1:100) = 0;
   R_isol(:,:,ii) = corr(B');
   absR_isol(:,:,ii) = corr(abs(B)');

end
toc

%% plot

figure;
colormap(broc.broc)




% == encoding profile

nexttile;
imagesc((condVec \ [Y1{1}, Y2{1}])', [-3,3])
axis('square')
% xlabel('conditions')
% ylabel('voxels')
title('R pure')


nexttile;
% y10 = Y1{1}; y10(:, 1:50) = 0;
% y20 = Y2{1}; y20(:, 1:50) = 0;
B = condVec \ [zeros(900,100), Y1{1}(:,51:end), Y2{1}(:,51:end)];
imagesc(B', [-3,3])
axis('square')
title('R pure isol')



nexttile;
imagesc((condVec \ [Ycomb{1}])', [-3,3])
axis('square')
% xlabel('conditions')
% ylabel('voxels')
title('R mixed')



nexttile;
B = (condVec \ [Ycomb{1}])';
B(1:100,:) = 0;
imagesc(B, [-3,3])
axis('square')
% xlabel('conditions')
% ylabel('voxels')
title('R mixed isol')




% == R

% nexttile;
% imagesc(mean(R_sep,3)./std(R_sep,[],3), [-1 1])
% axis('square')
% title('R pure')
% 
% 
% nexttile;
% imagesc(mean(R_isol_sep,3)./std(R_isol_sep,[],3), [-1 1])
% axis('square')
% title('R pure isol')
% 
% 
% 
% nexttile;
% imagesc(mean(R_comb,3)./std(R_comb,[],3), [-1 1])
% axis('square')
% title('R mixed')
% 
% 
% nexttile;
% imagesc(mean(R_isol,3)./std(R_isol,[],3), [-1 1])
% axis('square')
% title('R mixed isol')



nexttile;
imagesc(mean(R_sep,3), [-1 1])
axis('square')
title('R pure')


nexttile;
imagesc(mean(R_isol_sep,3), [-1 1])
axis('square')
title('R pure isol')



nexttile;
imagesc(mean(R_comb,3), [-1 1])
axis('square')
title('R mixed')


nexttile;
imagesc(mean(R_isol,3), [-1 1])
axis('square')
title('R mixed isol')







% == absR


% 
% nexttile;
% imagesc(mean(absR_sep,3)./std(absR_sep,[],3), [-20 20])
% axis('square')
% title('absR pure')
% 
% 
% 
% nexttile;
% imagesc(mean(absR_isol_sep,3)./std(absR_isol_sep,[],3), [-20 20])
% axis('square')
% title('R isol pure')
% 
% 
% 
% nexttile;
% imagesc(mean(absR_comb,3)./std(absR_comb,[],3), [-20 20])
% axis('square')
% title('absR mixed')
% 
% 
% nexttile;
% imagesc(mean(absR_isol,3)./std(absR_isol,[],3), [-20 20])
% axis('square')
% title('absR mixed isol')




nexttile;
imagesc(mean(absR_sep,3), [-1 1])
axis('square')
title('absR pure')



nexttile;
imagesc(mean(absR_isol_sep,3), [-1 1])
axis('square')
title('R isol pure')



nexttile;
imagesc(mean(absR_comb,3), [-1 1])
axis('square')
title('absR mixed')


nexttile;
imagesc(mean(absR_isol,3), [-1 1])
axis('square')
title('absR mixed isol')













%% RSA sim sample

bam = load('/Users/hr0283/Dropbox (Brown)/RDM_fmri_scripts/util/ScientificColourMaps7/bam/bam.mat')
center = @(x) x-nanmean(x);
vec = @(x) x(:);

noises = 10.^linspace(-1,1,10);
clear ega_r rsa_r ega_d rsa_d enc_d
tic
parfor nn = 1:length(noises)


    %% parameters
    nTrials = 150;
    D = struct;
    D.numPart = 1;
    D.numVox  = 250;
    numSim = 1000;
    signal = .1;
    noise  = noises(nn);

    %% model

    rel1a = 1;
    rel1b = 0;

    rel2a = 0;
    rel2b = 0;


    M  = struct;
    M.type       = 'feature';
    M.numGparams = 2;
    M.Ac(:,:,1) = [1 0 rel1a 0; 0 0 0 0; rel1b 0 1 0; 0 0 0 0];
    M.Ac(:,:,2) = [0 0 0 0; 0 1 0 rel2a;  0 0 0 0; 0 rel2b 0 1];
    M.theta0=ones(2,1);                        % Starting values: could be closer, but converges anyways


    %% design matrix
    xs = repmat([-2:2]', [nTrials/5, 1]);

    X1 = [xs(randperm(nTrials)), xs(randperm(nTrials))];
    X2 = [xs(randperm(nTrials)), xs(randperm(nTrials))];
    X = blkdiag(X1 - mean(X1), X2 - mean(X2));



    X1_rsa = [];
    X2_rsa = [];
    pred1_rsa = [];
    pred2_rsa = [];


    for ii = [-2, -1, 1, 2]
        for jj = [-2, -1, 1, 2]

            X1_rsa = [X1_rsa, (X1(:,1) == ii) & (X1(:,2) == jj)];
            X2_rsa = [X2_rsa, (X2(:,1) == ii) & (X2(:,2) == jj)];

            pred1_rsa = [pred1_rsa; ii];
            pred2_rsa = [pred2_rsa; jj];


        end
    end

    X_rsa = blkdiag(X1_rsa - mean(X1_rsa), X2_rsa - mean(X2_rsa));

    rsa_triu = triu(true(16),1);

    rsa_pred1 = abs(pred1_rsa - pred1_rsa');
    v_pred1 = rsa_pred1(rsa_triu);

    rsa_pred2 = abs(pred2_rsa - pred2_rsa');
    v_pred2 = rsa_pred2(rsa_triu);


    %% generate
    [Y, partVec, condVec] = pcm_generateData(M, M.theta0, D, numSim, signal, noise, 'design', X);



    %% get EGA similarity
    R = zeros(2,numSim);
    for ii = 1:length(Y)

        B = condVec \ Y{ii};
        r = corr(B(1:2,:)', B(3:4,:)', 'rows', 'pairwise');
        R(:,ii) = [r(1,1), r(2,2)];

    end

    %% get RSA similarity
    R_rsa = zeros(2, numSim);
    for ii = 1:length(Y)

        B_rsa = X_rsa \ Y{ii};
        r_rsa = corr(B_rsa(1:16,:)', B_rsa(17:32,:)', 'rows', 'pairwise');
        %         R_rsa(:,ii) = corr([v_pred1, v_pred2], 1-r_rsa(rsa_triu), 'rows', 'pairwise');
        b = center([v_pred1, v_pred2]) \ center(1-r_rsa(rsa_triu));
        %         b = ([v_pred1*0 + 1, v_pred1, v_pred2]) \ (1-r_rsa(rsa_triu));
        R_rsa(:,ii) = b(1:2);

    end


    %% get encoding
    R_enc = zeros(2, numSim);
    for ii = 1:length(Y)

        B_enc = condVec(1:150,1:2) \ Y{ii}(1:150,:);
        ypred1 = condVec(151:300,3) * B_enc(1,:);
        ypred2 = condVec(151:300,4) * B_enc(2,:);

        R_enc(:,ii) = center([vec(ypred1), vec(ypred2)]) \ center(vec(Y{ii}(151:300, :)));

    end



    %% save

    ega_r(:,nn) = mean(R,2);
    rsa_r(:,nn) = mean(R_rsa,2);

    ega_d(:,nn) = mean((R),2) ./ std((R),[],2);
    rsa_d(:,nn) = mean((R_rsa),2) ./ std((R_rsa),[],2);
    enc_d(:,nn) = mean((R_enc),2) ./ std((R_enc),[],2);


%     ega_d(:,nn) = corr(R', signal) 
%     rsa_d(:,nn) = corr(R_rsa', signal);
%     enc_d(:,nn) = corr(R_enc', signal);



    %% plot

    % nexttile; hold on;
    % ksdensity(R(1,:))
    % ksdensity(R_rsa(1,:))
    % xline(rel1)
    % legend({'EGA', 'RSA'})
    % xlim([0,1]);
    %
    % nexttile; hold on;
    % ksdensity(R(2,:))
    % ksdensity(R_rsa(2,:))
    % xline(0)
    % legend({'EGA', 'RSA'})
    % % xlim([-.25,.25]);


end
toc


%% plot
figure;

% nexttile; hold on;
% plot(log10(noises), ega_r(1,:)', '-r', 'LineWidth', 2)
% plot(log10(noises), rsa_r(1,:)', '-b', 'LineWidth', 2)
% title('mean R - signal')
% axis('square')
% set(gca, 'TickDir', 'out', 'LineWidth', 1)
% legend({'EGA', 'RSA'})
%
%
% nexttile; hold on;
% plot(log10(noises), ega_r(2,:)', '-r', 'LineWidth', 2)
% plot(log10(noises), rsa_r(2,:)', '-b', 'LineWidth', 2)
% title('mean R - noise')
% axis('square')
% set(gca, 'TickDir', 'out', 'LineWidth', 1)



nexttile; hold on;
plot(log10(noises), ega_d(1,:)', '-r', 'LineWidth', 2)
plot(log10(noises), rsa_d(1,:)', '-b', 'LineWidth', 2)
plot(log10(noises), enc_d(1,:)', '-g', 'LineWidth', 2)
title('cohen d - signal')
axis('square')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
legend({'EGA', 'RSA', 'Encoding'})


nexttile; hold on;
plot(log10(noises), ega_d(2,:)', '-r', 'LineWidth', 2)
plot(log10(noises), rsa_d(2,:)', '-b', 'LineWidth', 2)
plot(log10(noises), enc_d(2,:)', '-g', 'LineWidth', 2)
title('cohen d - noise')
axis('square')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
legend({'EGA', 'RSA', 'Encoding'})






