
% Generate signal
srate = 1000;

frex   = [ 3 10 5 15 35];
amplit = [ 5 15 10 5 7];

time = -1:1/srate:1;

signal = zeros(1, length(time));
for fi=1:length(frex)
    signal = signal + amplit(fi) * sin(2*pi*frex(fi)*time);
end

% Slow FT in a loop
N  = length(signal);
fourierTime = (0:N-1)/N;
nyquist = srate/2;

fourierCoefs = zeros(size(signal));
frequencies = linspace(0, nyquist, floor(N/2) + 1);

for fi=1:N
    fourierSine = exp(-1i * 2 * pi * fourierTime * (fi - 1));

    fourierCoefs(fi) = dot(fourierSine, signal);
end

% normalize
fourierCoefs = fourierCoefs / N;

% plot
figure(1), clf
subplot(221)
plot(real(exp(-2*pi*1i*(10).*fourierTime)))
xlabel('time (a.u.)', ylabel('Amplitude'))
title('One sine wave from the FT (real part)')

subplot(222)
plot(signal)
title('Data')

subplot(212)
plot(frequencies, abs(fourierCoefs(1:length(frequencies)))*2, '*-')
xlabel('Frequency (Hz)'), ylabel('Amplitude (\muV)')
title('Amplitude spectrum derived from discrete Fourier transform')

% Fast FT
fourierCoefsF = fft(signal) / N;

subplot(212), hold on
plot(frequencies, abs(fourierCoefsF(1:length(frequencies)))*2, 'ro')
set(gca, 'xlim', [0 40])
