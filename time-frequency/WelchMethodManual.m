load EEGrestingState.mat
N = length(eegdata);

% define the time vector
timevec = (0:N-1)/srate;

% plot the data
figure(1), clf
plot(timevec, eegdata, 'k')
xlabel('Time (seconds)'), ylabel('Voltage (\muV)')

% one big FT
eegpow = abs( fft(eegdata)/N ).^2;
hz = linspace(0, srate/2, floor(N/2)+1);

% manual Welch method
winLength = 1*srate;
nSkip = round(srate/2);

% window onset times
winonsets = 1:nSkip:N-winLength;
hzW = linspace(0, srate/2, floor(winLength/2)+1);

% Hann win
hannw = .5 - cos(2*pi*linspace(0, 1, winLength))./2;

% initialize power matrix (windows x frequencies)
eegpowW = zeros(1, length(hzW));

% loop over frequencies
for wi=1:length(winonsets)
    % get a chunck of data
    datachunk = eegdata(winonsets(wi):winonsets(wi)+winLength-1);
    datachunk = datachunk .* hannw;

    tmppow = abs(fft(datachunk)/winLength).^2;
    eegpowW = eegpowW + tmppow(1:length(hzW));
end

% plotting
figure(2), clf
plot(hz, eegpow(1:length(hz)), 'k', 'linew', 2)
hold on
plot(hzW, eegpowW/2000, 'r', 'linew', 2)
set(gca, 'xlim', [0 40])
xlabel('Frequency (Hz)')
legend({'FFT';'Welch'})
title('Power spectrum (Welch method)')


%% MATLAB pwelch

subplot(212)

% create Hann window
winsize = 2*srate; % 2-second window
hannw = .5 - cos(2*pi*linspace(0,1,winsize))./2;

% number of FFT points (frequency resolution)
nfft = srate*100;

pwelch(eegdata,hannw,round(winsize/4),nfft,srate);
set(gca,'xlim',[0 40])


%   Welch's method on v1 laminar data


clear
load v1_laminar.mat
csd = double(csd);

% specify a channel for the analyses
chan2use = 7;


% create Hann window
hannw = .5 - cos(2*pi*linspace(0,1,size(csd,2)))./2;

% Welch's method using MATLAB pwelch
[pxx,hz] = pwelch(squeeze(csd(chan2use,:,:)),hannw,round(size(csd,2)/10),1000,srate);

figure(28), clf
subplot(211)
plot(timevec,mean(csd(chan2use,:,:),3),'linew',2)
set(gca,'xlim',timevec([1 end]))
xlabel('Time (s)'), ylabel('Voltage (\muV)')

subplot(212)
plot(hz,mean(pxx, 2),'linew',2)
set(gca,'xlim',[0 140])
xlabel('Frequency (Hz)')
ylabel('Power (\muV^2)')
