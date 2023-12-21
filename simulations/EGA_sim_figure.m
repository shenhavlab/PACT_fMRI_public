%% generate figure

close all;

nnoise = logspace(-1,1, 10);
nsim = 1000;
nsim_null = 1000;
vec = @(x) x(:);

% original data
nFeat = 250;
nTrials = 1000;

B0_train = randn(nFeat,1);
X_train = randn(nTrials,1);
X_test = randn(nTrials,1);


figure; hold on;
tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact')
cols = colormap('winter');
cidx = round(linspace(1,size(cols,1)*.8, length(nnoise)));





% ==== sensitivity
clear rmsCB rmsCY rmsTB

for ss = 1:length(nnoise)

    [corr_B, corr_Y, true_B] = deal(nan(nsim,1));

 

    parfor ii = 1:nsim

        B0_test  = B0_train + randn(nFeat,1);

        Y0_train = X_train*B0_train' + randn(nTrials,nFeat)*nnoise(ss); % match noise in ground truth and simulation
        Y0_test  = X_test*B0_test' + randn(nTrials,nFeat)*nnoise(ss); % match noise in ground truth and simulation

        B_train = X_train\Y0_train;
        B_test = X_test\Y0_test;

        Y_test = X_test*B_train;


        corr_B(ii) = corr(B_train', B_test');
        corr_Y(ii) = corr(vec(Y0_test), vec(Y_test));
        true_B(ii) = corr(B0_train, B0_test);
%         true_B(ii) = corr(vec(Y0_train), vec(Y0_test));

    end

    % EGO a
    nexttile(2); hold on;
    plot(corr_B, corr_Y, 'o', 'Color', cols(cidx(ss),:), 'LineWidth', 1.5, 'MarkerSize', 10);
    xlabel('beta corr')
    ylabel('y corr')
    plot([0,1], [0,1], '-k', 'LineWidth',1)
    set(gca, 'TickDir', 'out', 'LineWidth', 1)
    axis('square')

    nexttile(1); hold on;
    plot(true_B, corr_B, 'o', 'Color', cols(cidx(ss),:), 'LineWidth', 1.5, 'MarkerSize', 10);
    title(sprintf('B noise ~ B corr'))
    xlabel('beta noise')
    ylabel('beta corr')
    plot([0,1], [0,1], '-k', 'LineWidth',1)
    set(gca, 'TickDir', 'out', 'LineWidth', 1)
    axis('square')

    rmsCB(ss) = corr(true_B, corr_B);
    rmsCY(ss) = corr(true_B, corr_Y);

end


% ==== unbiased

[corr_B, corr_Y, s_B, s_Y, true_B] = deal(nan(nsim,1));

parfor ii = 1:nsim_null

    B0_test  = randn(nFeat,1);

    Y0_train = X_train*B0_train' + randn(nTrials,nFeat)*nnoise(1); % match noise in ground truth and simulation
    Y0_test  = X_test*B0_test' + randn(nTrials,nFeat)*nnoise(1); % match noise in ground truth and simulation

    B_train = X_train\Y0_train;
    B_test = X_test\Y0_test;

    Y_test = X_test*B_train;


    corr_B(ii) = corr(B_train', B_test');
    corr_Y(ii) = corr(vec(Y0_test), vec(Y_test));
    true_B(ii) = corr(B0_train, B0_test);

end

nexttile(4); hold on;
[f,x]=ksdensity(corr_B);
plot(x,f, '-k', 'LineWidth', 2)
[f,x]=ksdensity(corr_Y);
plot(x,f, '-r', 'LineWidth', 2)
xline(0)
title('unbiased')
legend({'EGA', 'encoding'})
set(gca, 'TickDir', 'out', 'LineWidth', 1)
axis('square')
xlabel('encoding strength under null')
ylabel('pdf')



% ==== comparision
nexttile(3); hold on;
plot(log10(nnoise), rmsCB, '-k', 'LineWidth', 2)
plot(log10(nnoise), rmsCY, '-r', 'LineWidth', 2)
legend({'EGA', 'encoding'})
title('high r')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
axis('square')
xlabel('log residual variance')
ylabel('correlation with ground-truth pattern similarity')


