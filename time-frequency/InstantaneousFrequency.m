% Instantaneous Frequency - a way to describe how the frequency of the signal changes over time.

load EEGrestingState.mat
time = (0:length(eegdata)-1)/srate;

figure(1), clf;
pwelch(eegdata, srate, srate/2, srate*2, srate);
set(gca, 'xlim', [0 60]);

% apply a narrowband filter at 10hz
nyquist = srate/2;
frange = [8 12];
order = round(20*srate/frange(1));

filtkern = fir1(order, frange/nyquist);

% compute the power spectrum of the filter kernel
filtpow = abs(fft(filtkern)).^2;
hz = linspace(0, nyquist, floor(length(filtkern)/2)+1);
filtpow = filtpow(1:length(hz));

% visualize the filter kernel
figure(2), clf;
subplot(121);
plot(filtkern, 'linew', 2);
xlabel('Time Points');
title('Filter kernel(fir1)');
axis square;

% plot amplitude spectrum of the filter kernel
subplot(122), hold on;
plot(hz, filtpow, 'ks-', 'linew', 2, 'markerfacecolor', 'w');
plot([0 frange(1) frange frange(2) nyquist], [0 0 1 1 0 0], 'ro-', 'linew', 2, 'markerfacecolor', 'w');

% dotted line corresponding to the lower edge of the filter cut-off
plot([1 1]*frange(1), get(gca, 'ylim'), 'k:');

% make the plot look nicer
set(gca, 'xlim', [0 frange(1)*4]);
xlabel('Frequency (Hz)'), ylabel('Filter gain');
legend({'Actual'; 'Ideal'});
title('Frequency response of filter (fir1)');
axis square;

% IF on filtered data
% apply 10hz filter
feegdata = filtfilt(filtkern, 1, double(eegdata));

% plot data for comparison
figure(3), clf;
plot(time, eegdata, time, feegdata, 'linew', 2);
xlabel('Time (s)'), ylabel('Amplitude');

% the actual instantaneous frequency part
angels = angle(hilbert(feegdata));
instalpha = diff(unwrap(angels)) / (2*pi/srate);

figure(4), clf;
plot(time(1:end-1), instalpha, 'ks-', 'markerfacecolor', 'w');
xlabel('Time (s)'), ylabel('Frequency (Hz)');

% apply median filter to supra-threshold points

% convert to z and show histogram
instz = (instalpha - mean(instalpha)) / std(instalpha);

figure(5), clf;
hist(instz, 200);
xlabel('I.F. (z)'), ylabel('Count');
title('Distribution of z-norm. IF');

% identify supra-threshold data points
tofilter = find(abs(instz)>2);

% median filter
instalphaFilt = instalpha;
k = round(1000*50/srate);
for i=1:length(tofilter)
    indices = max(1, tofilter(i)-k):min(pnts, tofilter(i)+k);
    instalphaFilt(tofilter(i)) = median(instalpha(indices));ßß
end

figure(6), clf, hold on
plot(time(1:end-1),instalpha,'k','linew',2)
plot(time(1:end-1),instalphaFilt,'r--','linew',2)
xlabel('Time (s)'), ylabel('Frequency (Hz)')
legend({'Original';'Filtered'})
set(gca,'ylim',[5 15])
