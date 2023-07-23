% reconstruct a signal from frequency domain
% define sampling rate and time vector
srate = 1000;
time = -1:1/srate:1;

% frequencies of the original signal
frex = [ 3 10 5 15 35 ];

signal = zeros(1, length(time));
for fi=1:length(frex)
    signal = signal + fi * sin(2 * pi * frex(fi) * time);
end

N = length(signal);
fourierTime = (0:N-1)/N;

reconSignal = zeros(size(signal));
fourierCoefs = fft(signal)/N;

for fi=1:N
    fourierSine = exp(1i * 2 * pi * fourierTime * (fi - 1));

    reconSignal = reconSignal + (fourierSine * fourierCoefs(fi));
end

figure (1), clf
plot(time, real(reconSignal), 'b-');
hold on
plot(time(1:10:end), signal(1:10:end), 'ro');
legend({'Reconstructed', 'Original'})
