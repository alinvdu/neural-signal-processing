load ../data/sampleEEGdata.mat

% print some details about the data
disp("Name of the dataset: " + EEG.setname);
chanLocsName = EEG.chanlocs(1).labels;

for k = 2:numel(EEG.chanlocs)
    chanLocsName = chanLocsName + ","  + EEG.chanlocs(k).labels;
end

disp("Number of channels: " + size(EEG.chanlocs, 2));
disp("Channels locations: " + chanLocsName);

% Time frequency analysis on multiple trials - 1 channel

% pick some frequency parameters
min_freq = 4;
max_freq = 30;
num_freq = 40;
frex = linspace(min_freq, max_freq, num_freq);

% pick a channel to compute
channel_to_use = 'o2';

% wavelets parameters
range_cycles = [4 10];
s = logspace(log10(range_cycles(1)), log10(range_cycles(2)), num_freq) ./ (2 * pi * frex); 
wave_time = -2:1/EEG.srate:2;
half_wave = (length(wave_time) - 1) / 2;

% FFT parameters
nWave = length(wave_time);
nData = EEG.pnts * EEG.trials;
nConv = nWave + nData - 1;

% compute FFT of all trials concatenated
all_data = reshape(EEG.data(strcmpi(channel_to_use, {EEG.chanlocs.labels}), :, :), 1, []);
dataX = fft(all_data, nConv);

% initialize output time frequency data
tf = zeros(num_freq, EEG.pnts);

% convolution over frequencies
for fi=1:length(frex)
    % create a wavelet and get its FFT
    wavelet = exp(2*1i*pi*frex(fi).*wave_time) .* exp(-wave_time.^2./(2*s(fi)^2));
    waveletX = fft(wavelet, nConv);

    % scale
    waveletX = waveletX ./ max(waveletX);

    % convolution
    as = ifft(waveletX .* dataX);

    % cut wings
    as = as(half_wave+1:end-half_wave);

    % reshape back to time by trial
    as = reshape(as, EEG.pnts, EEG.trials);

    % compute power and average over trials
    tf(fi, :) = mean(abs(as).^2, 2);
end

% plot
figure(1), clf;
contourf(EEG.times, frex, tf, 40, "linecolor", "none");
set(gca, "clim", [1 4], "xlim", [-500 1300]);
colormap hot

% figure meta data
title("Time & frequency plot of channel " + channel_to_use + " across all trials.");
xlabel("Time (ms)");
ylabel("Frequency (Hz)");