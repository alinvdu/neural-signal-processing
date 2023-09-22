% Inter Trial Phase Clustering - quantifying the consistency of phase of neural signals across trials.
% Sunchronization check of neural activity between trials.

load sampleEEGdata.mat

% wavelet parameters
num_frex = 40;
min_freq = 2;
max_freq = 30;

channel2use = 'pz';

% set range for variable number of wavelet cycles
range_cycles = [ 3 10 ];

% other params
frex = logspace(log10(min_freq), log10(max_freq), num_frex);
nCycs = logspace(log10(range_cycles(1)), log10(range_cycles(end)), num_frex);
time = -2:1/EEG.srate:2;

% FFT params
nWave = length(time);
nDdata = EEG.pnts + EEG.trials;
nConv = nWave + nData - 1;

% FFT of data - it doesn't change so we can compute it outside the loop
dataX = fft(reshape(EEG.data(strcmpi(channel2use, {EEG.chanlocs.labels}), :, :), 1, nData), nConv);

% initialize output time-frequency data
tf = zeros(num_frex, EEG.pnts);

% loop over frequencies
for fi=1:num_frex
    % create wavelet and get its FFT
    s = nCycs(fi) / (2*pi*frex(fi));
    wavelet = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2 ./ (2*s^2));
    waveletX = fft(wavelet, nConv);

    % run convolution
    as = ifft(waveletX .* dataX, nConv);
    as = as(half_wave+1:end-half_wave);

    % reshape back
    as = reshape(as, EEG.pnts, EEG.trials);

    % compute ITPC
    tf(fi, :) = abs(mean(exp(1i * angle(as)), 2)); % phase angle checkc
end

% plot visualization
figure(1), clf
contourf(EEG.times, frex, tf, 40, 'linecolor', 'none');
set(gca, 'clim', [0 .6], 'ydir', 'normal', 'xlim', [-300 1000]);
title('ITPC');
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
