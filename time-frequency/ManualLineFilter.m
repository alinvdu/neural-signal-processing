% Generate 1/f noise with 50 hz line noise.
srate = 1234;
npnts = srate * 3;
time = (0:npnts-1)/srate;

ed = 50;
as = rand(1, npnts) .* exp(-(0:npnts-1)/ed);
fc = as .* exp(1i * 2 * pi * rand(size(as)));

signal = real(ifft(fc)) * npnts;
signal = signal + sin(2 * pi * 50 * time);

hz = linspace(0, srate/2, floor(npnts/2)+1);

% plot original signal and power spectrum
signalX = fft(signal);
figure(1), clf
subplot(211)
plot(time,signal,'k')
xlabel('Time (s)'), ylabel('Activity')
title('Time domain')

subplot(212)
plot(hz,(2*abs(signalX(1:length(hz)))/npnts).^2,'k')
set(gca,'xlim',[0 80])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Frequency domain')

% zero out the 50 hz component
hz50idx = dsearchn(hz',50);
signalXf = signalX;

signalXf(hz50idx) = 0;
signalXf(end-hz50idx+2) = 0; % zero out the negative frequency as well

signalff = ifft(signalXf);
signalXff = fft(signalff);

figure(2), clf
subplot(211), hold on
plot(time,signal,'k')
plot(time,signalff,'r')
xlabel('Time (s)'), ylabel('Activity')
legend({'Original';'Filtered'})
title('Time domain')

subplot(212), hold on
plot(hz,2*abs(signalX(1:length(hz)))/npnts,'k')
plot(hz,2*abs(signalXff(1:length(hz)))/npnts,'ro-','markerfacecolor','r')
set(gca,'xlim',[0 80])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Frequency domain')
