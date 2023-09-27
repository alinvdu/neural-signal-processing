% Phase synchronization in voltage and laplacian eeg data
load ../data/sampleEEGdata.mat

chan1 = 'FCz';
chan2 = 'POz';

% frequency params
min_freq = 2;
max_freq = 40;
num_frex = 34;

fwhm = linspace(.3, .1, num_frex);

% wavelet params
frex = logspace(log10(min_freq), log10(max_freq), num_frex);
time = -2:1/EEG.srate:2;
half_wave = (length(time)-1)/2;

% FFT params
nWave = length(time);
nData = EEG.pnts * EEG.trials;
nConv = nWave + nData - 1;

% FFT of data
data1X = fft(reshape(EEG.data(strcmpi(chan1, {EEG.chanlocs.labels}), :, :), 1, nData), nConv);
data2X = fft(reshape(EEG.data(strcmpi(chan2, {EEG.chanlocs.labels}), :, :), 1, nData), nConv);

% Inter Site Phase Clustering & Phase Index
[ispc, pli] = deal(zeros(num_frex, EEG.pnts));

for fi=1:num_frex
    % create wavelet and get its FFT
    wavelet = exp(2*1i*pi*frex(fi).*time) .* exp(-4*log(2)*time.^2 / fwhm(fi).^2);
    waveletX = fft(wavelet, nConv);
    waveletX = waveletX ./ max(waveletX);

    % convolution for chan1
    as1 = ifft(waveletX.*data1X, nConv);
    as1 = as1(half_wave+1:end-half_wave);
    as1 = reshape(as1, EEG.pnts, EEG.trials);

    % convolution for chan2
    as2 = ifft(waveletX.*data2X, nConv);
    as2 = as2(half_wave+1:end-half_wave);
    as2 = reshape(as2, EEG.pnts, EEG.trials);

    % collect phase angle differences
    cdd = exp(1i*(angle(as1) - angle(as2)));

    % compute ISPC and PLI (and average over trials)
    ispc(fi, :) = abs(mean(cdd, 2));
    pli(fi, :) = abs(mean(sign(imag(cdd)), 2));
end

% plotting
figure(1), clf;

subplot(221);
contourf(EEG.times, frex, ispc, 40, 'linecolor', 'none');
set(gca, 'xlim', [-300 1200], 'clim', [0 .4]);
colormap hot; colorbar;
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title('ISCP, Voltage');

subplot(222);
contourf(EEG.times, frex, pli, 40, 'linecolor', 'none');
set(gca, 'xlim', [-300 1200], 'clim', [0 .4]);
colormap hot; colorbar;
title('PLI, Voltage');

% repeat procedure for Laplacian
% compute laplacian - store in another variable
EEG.lap = laplacian_perrinX(EEG.data, [EEG.chanlocs.X], [EEG.chanlocs.Y], [EEG.chanlocs.Z]);

% FFT of data
data1X = fft(reshape(EEG.lap(strcmpi(chan1, {EEG.chanlocs.labels}), :, :), 1, nData), nConv);
data2X = fft(reshape(EEG.lap(strcmpi(chan2, {EEG.chanlocs.labels}), :, :), 1, nData), nConv);

% Inter Site Phase Clustering & Phase Index
[ispc, pli] = deal(zeros(num_frex, EEG.pnts));

for fi=1:num_frex
    % create wavelet and get its FFT
    wavelet = exp(2*1i*pi*frex(fi).*time) .* exp(-4*log(2)*time.^2 / fwhm(fi).^2);
    waveletX = fft(wavelet, nConv);
    waveletX = waveletX ./ max(waveletX);

    % convolution for chan1
    as1 = ifft(waveletX.*data1X, nConv);
    as1 = as1(half_wave+1:end-half_wave);
    as1 = reshape(as1, EEG.pnts, EEG.trials);

    % convolution for chan2
    as2 = ifft(waveletX.*data2X, nConv);
    as2 = as2(half_wave+1:end-half_wave);
    as2 = reshape(as2, EEG.pnts, EEG.trials);

    % collect phase angle differences
    cdd = exp(1i*(angle(as1) - angle(as2)));

    % compute ISPC and PLI (and average over trials)
    ispc(fi, :) = abs(mean(cdd, 2));
    pli(fi, :) = abs(mean(sign(imag(cdd)), 2));
end

subplot(223);
contourf(EEG.times, frex, ispc, 40, 'linecolor', 'none');
set(gca, 'xlim', [-300 1200], 'clim', [0 .4]);
colormap hot; colorbar;
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title('ISCP, Laplacian');

subplot(224);
contourf(EEG.times, frex, pli, 40, 'linecolor', 'none');
set(gca, 'xlim', [-300 1200], 'clim', [0 .4]);
colormap hot; colorbar;
title('PLI, Laplacian');
