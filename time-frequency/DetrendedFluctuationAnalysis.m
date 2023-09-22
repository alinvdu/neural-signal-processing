load dfa_data.mat

% Steps
% 1. Convert to mean-centered cum sum.
% 2. Define log-spaced scales.
% 3. Split signal into epochs, detrend, compute RMS.
% 4. Repeat (3) for each scale.
% 5. Compute linear regression between log scales andd log-RMSes.

xfilt = abs(hilbert(filterFGx(x, srate, 10, 5)));

% create ddata with DFA = .5 -> random noise has DFA .5
N = length(x);
randnoise = randn(N, 1);

% parameters
nScales = 20;
ranges = round(N*[.01 .2]);
scales = ceil(logspace(log10(ranges(1)), log10(ranges(2)), nScales));
rmses = zeros(2, nScales);

% plot the two signals
figure(1), clf
subplot(221)
plot(timevec, randnoise)
title('Signal 1: white noise')
xlabel('Time (seconds)')

subplot(222)
plot(timevec, xfilt)
title('Signal 2: real data')
xlabel('Time (seconds)')

% step 1:
randnoise = cumsum(randnoise - mean(randnoise));
x4dfa = cumsum(xfilt - mean(xfilt));

% show time series for comparison
subplot(223)
plot(timevec, randnoise)
title('Integrated noise')

subplot(224)
plot(timevec, x4dfa)
title('Integrated signal')

% step 3 and 4
for scalei=1:nScales
    % number of epochs for this scale
    n = floor(N/scales(scalei));

    % compute RMS for random noise
    epochs = reshape(randnoise(1:n*scales(scalei)), scales(scalei), n);
    depochs = detrend(epochs);

    % rms computation
    rmses(1, scalei) = mean(sqrt(mean(depochs.^2, 1)));

    % repeat for the real data
    epochs = reshape(x4dfa(1:n*scales(scalei)), scales(scalei), n);
    depochs = detrend(epochs);

    % rms computation
    rmses(2, scalei) = mean(sqrt(mean(depochs.^2, 1)));
end

% fit a linear model to quantify
A = [ones(nScales, 1) log10(scales)']; % linear model
dfa1 = (A' * A) \ (A'*log10(rmses(1, :))'); % fit to noise
dfa2 = (A' * A) \ (A'*log10(rmses(2, :))'); % fit to signal

% plot the linear fit (in log-log space)
figure(12), clf, hold on

% plot results for white noise
plot(log10(scales), log10(rmses(1, :)), 'rs', 'linew', 2, 'markerfacecolor', 'w', 'markersize', 10)
plot(log10(scales), dfa1(1) + dfa1(2) * log10(scales), 'r--', 'linew', 2)

% plot results for the real signal
plot(log10(scales), log10(rmses(2, :)), 'bs', 'linew', 2, 'markerfacecolor', 'w', 'markersize', 10)
plot(log10(scales), dfa2(1) + dfa2(2) * log10(scales), 'b--', 'linew', 2)

legend({'Data (noise)';[ 'Fit (DFA=' num2str(round(dfa1(2),3)) ')' ]; ...
        'Data (signal)';[ 'Fit (DFA=' num2str(round(dfa2(2),3)) ')' ] })
xlabel('Data scale (log)'), ylabel('RMS (log)')
title('Comparison of Hurst exponent for different noises')
axis square
