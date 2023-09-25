% Measure connection in terms of phase between two electrodes on trial 1
load ../data/v1_laminar

chan1idx = 1;
chan2idx = 8;

% complex morlet wavelet
cent_freq = 8;
time = -1.5:1/srate:1.5;
s = 8/(2*pi*cent_freq);
wavelet = exp(2*1i*pi*cent_freq.*time) .* exp(-time.^2/(2*s^2));
half_wavN = (length(time)-1)/2;

% FFT params
nWave = length(time);
nData = size(csd, 2);
nConv = nWave + nData - 1;

%% Perform FFT on wavelet
waveletX = fft(wavelet, nConv);
waveletX = waveletX ./ max(waveletX);

% initialize outtput time-frequency data
phase_data = zeros(2, nData);
real_data = zeros(2, nData);

% analytic signal of channel 1
dataX = fft(csd(chan1idx, :, 1), nConv);
as = ifft(waveletX .* dataX, nConv);
as = as(half_wavN+1:end-half_wavN);

phase_data(1, :) = angle(as);
real_data(1, :) = real(as);

% analytic signal of channel 8
dataX = fft(csd(chan2idx, :, 1), nConv);
as = ifft(waveletX .* dataX, nConv);
as = as(half_wavN+1:end-half_wavN);

% collect real and phase data
phase_data(2, :) = angle(as);
real_data(2, :) = real(as);

% quantify phase synch between the two channels
phase_synchronization = abs(mean(exp(1i*diff(phase_data))));

disp(['Synchronization between ' num2str(chan1idx) ' and ' num2str(chan2idx) ' is ' num2str(phase_synchronization) '!'])
figure(1), clf;
h = polar([0 angle(mean(exp(1i*diff(phase_data))))], [0 phase_synchronization]);
set(h, 'linewidth', 6, 'color', 'g');
