% Plot different baseline normalization

clear
load sampleEEGdata.mat

baseline_windows = [ -500 200;
                     -100   0;
                        0 300;
                     -800   0;];

% convert baseline time to indices
baseidx = reshape(dsearchn(EEG.times', baseline_windows(:)), [], 2);

% setup wavelet params
min_freq = 2;
max_freq = 30;
num_frex = 40;
frex = linspace(min_freq, max_freq, num_frex);

channel2use = 'o1';
range_cycles = [ 4 10 ];

s = logspace(log10(range_cycles(1)), log10(range_cycles(end)), num_frex) ./ (2*pi*frex);
wavtime = -2:1/EEG.srate:2;
half_wave = (length(wavtime)-1)/2;

% FFT Params
nWave = length(wavtime);
nData = EEG.pnts * EEG.trials;
nConv = nWave + nData - 1;

alldata = reshape(EEG.data(strcmpi(channel2use, {EEG.chanlocs.labels}), :, :), 1, []);
dataX = fft(alldata, nConv);

% initialize time-frequency for TF
tf = zeros(size(baseidx, 1), length(frex), EEG.pnts);

% convolution
for fi=1:length(frex)
    % create wavelet and get its FFT
    wavelet = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
    waveletX = fft(wavelet, nConv);
    waveletX = waveletX ./ max(waveletX);

    as = ifft(waveletX .* dataX);
    as = as(half_wave+1:end-half_wave);

    as = reshape(as, EEG.pnts, EEG.trials);

    tf(4, fi, :) = mean(abs(as).^2, 2);
end

% db normalization
clim = [-3 3];

% create a new matrix for percent change
tfpct = zeros(size(tf));

for basei=1:size(tf, 1)
    activity = squeeze(tf(4, :, :));
    baseline = mean(tf(4, :, baseidx(basei, 1):baseidx(basei, 2)), 3);

    % decibel
    tf(basei, :, :) = 10 * log10(activity ./ repmat(baseline', [1 EEG.pnts]));
end

figure(1), clf
for basei=1:size(baseline_windows,1)
    
    subplot(2,2,basei)
    
    contourf(EEG.times,frex,squeeze(tf(basei,:,:)),40,'linecolor','none')
    set(gca,'clim',clim,'ydir','normal','xlim',[-300 1000])
    title([ 'DB baseline of ' num2str(baseline_windows(basei,1)) ' to ' num2str(baseline_windows(basei,2)) ' ms' ])
end

xlabel('Time (ms)'), ylabel('Frequency (Hz)')
