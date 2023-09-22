load sampleEEGdata.mat

% extract reaction time from each trial in the EEGLAB data format.
rts = zeros(size(EEG.epoch));

for ei=1:EEG.trials
    % find the index corresponding to time=0, i.e., trial onset
    [~, zeroloc] = min(abs(cell2mat(EEG.epoch(ei).eventlatency)));

    % reaction time is the event after the trial onset
    rts(ei) = EEG.epoch(ei).eventlatency{zeroloc+1};
end

figure(1), clf
plot(rts, 'ks-', 'markerfacecolor', 'w', 'markersize', 12);
xlabel('Trial'), ylabel('Reaction time (ms)');

% Create the design matrix -> two regressors: intercept and RTs
X = [ ones(EEG.trials, 1) rts'];

freqrange = [2 25];
numfrex = 30;

% setup convolution parameters
wavtime = -2:1/EEG.srate:2;
frex = linspace(freqrange(1), freqrange(2), numfrex);
nData = EEG.pnts * EEG.trials;
nKern = length(wavtime);
nConv = nData + nKern - 1;
halfwav = (length(wavtime)-1)/2;
nCyc = logspace(log10(4), log10(12), numfrex);

% initialize time-frequency matrix
tf3d = zeros(numfrex, EEG.pnts, EEG.trials);

% compute Fourier coeff of EEG data - this can be computed outside of freq loop
eegX = fft(reshape(EEG.data(47, :, :), 1, []), nConv);

for fi=1:numfrex
    % create the wavelet
    s = nCyc(fi) / (2*pi*frex(fi));
    cmw = exp(2*1i*pi*frex(fi).*wavtime) .* exp((-wavtime.^2) / (2*s.^2));
    cmwX = fft(cmw, nConv);
    cmwX = cmwX ./ max(cmwX);

    as = ifft(eegX .* cmwX);

    % cut wings
    as = as(halfwav+1:end-halfwav);
    as = reshape(as, EEG.pnts, EEG.trials);

    % extract power
    tf3d(fi, :, :) = abs(as).^2;
end

% plot the time frequency plots of some trials
for i=1:3
    subplot(2, 3, i);
    imagesc(EEG.times, frex, squeeze(tf3d(:, :, i)));
    axis square, axis xy;
    set(gca, 'clim', [0 10], 'xlim', [-200 1200])
    xlabel('Time (ms)'), ylabel('Frequency');
    title(['Trial' num2str(i)]);
end

% show the trial-average map
subplot(212);
imagesc(EEG.times, frex, squeeze(mean(tf3d, 3)));
axis square, axis xy;
set(gca, 'clim', [0 5], 'xlim', [-200 1200]);
xlabel('Time (ms)'), ylabel('Frequency');
title('All trials');

% setup the regression model
% reshape the 3D matrix to 2D so we concat some parts together as shortcut
tf2d = reshape(tf3d, numfrex * EEG.pnts, EEG.trials)';

% fit the model on the 2D matrix
b = (X'*X)\X'*tf2d; % linear regression fit with OLS (Ordinary Least Squares)

% reshape b into a time by frequency matrix
betamat = reshape(b(2, :), length(frex), EEG.pnts);

% show the design and data matrices (you won't understand anything from them)
figure(2), clf

ax1_h = axes;
set(ax1_h,'Position',[.05 .1 .1 .8])
imagesc(X);
set(ax1_h,'xtick',1:2,'xticklabel',{'Int';'RTs'},'ydir','norm')
ylabel('Trials');
title('Design matrix');

ax2_h = axes;
set(ax2_h,'Position',[.25 .1 .7 .8])
imagesc(tf2d);
set(ax2_h,'ydir','norm','clim',[0 20])
ylabel('Trials');
xlabel('Time Frequency');
title('Data Matrix');

colormap gray;

% results - remember each point color intensity relates to the regression model, more intensity more correlation

figure(3), clf

contourf(EEG.times, frex, betamat, 40, 'linecolor', 'none');
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
set(gca, 'xlim', [-200 1200], 'clim', [-.012 .012]);
title('Regression against RT over trials');
