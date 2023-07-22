% Topological plot of alpha power avg over the scalp
clear
load restingstate64chans.mat

% convert to double precision
EEG.data = double(EEG.data);

% extract power
chanpowr = (2 * abs(fft(EEG.data(), [], 2)/EEG.pnts)).^2;

% average over trials
chanpowr = mean(chanpowr, 3);

% vector of frequencies
hz = linspace(0, EEG.srate/2, floor(EEG.pnts/2)+1);

% plot power spectrum of all channels
figure(1), clf
plot(hz, chanpowr(:, 1:length(hz)), 'linew', 2);
xlabel('Frequency (Hz)'), ylabel('Power (\muV^2)');
set(gca, 'xlim', [0 30], 'ylim', [0 50]);

% extract alpha power
alphabounds = [7 13];

% convert to indices
freqidx = dsearchn(hz', alphabounds');

% extract average alpha power
alphapower = mean(chanpowr(:, freqidx(1):freqidx(2)), 2);

% and plot
figure(2), clf
topoplotIndie(alphapower, EEG.chanlocs, 'numcontour', 0);
set(gca, 'clim', [0 6]);
colormap hot
