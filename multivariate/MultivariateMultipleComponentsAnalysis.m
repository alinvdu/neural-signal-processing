% Multivariate Multiple component analysis.
% Analyse first 2 components - compute power spectrum.
% Analyse all 8 components from GED and perform correlation analysis.

load('../data/restingstate.mat');
plot_simEEG(EEG, 31, 10);
EEG.fdata = filterFGx(EEG.data, EEG.srate, 11, 7);

% generate S and R covariate matrices
tidx = dsearchn(EEG.times', [0 2000]');
[S, R] = deal(zeros(EEG.nbchan));

for triali=1:EEG.trials
    % R matrix
    tmp = EEG.data(:, tidx(1):tidx(2), triali);
    tmp = bsxfun(@minus, tmp, mean(tmp, 2));
    R = R + tmp*tmp'/diff(tidx);

    % S matrix
    tmp = EEG.fdata(:, tidx(1):tidx(2), triali);
    tmp = bsxfun(@minus, tmp, mean(tmp, 2));
    S = S + tmp*tmp'/diff(tidx);
end

S = S/EEG.trials;
R = R/EEG.trials;

clim = [-1 1]*max(abs([S(:); R(:)]))*.7;

figure(1), clf; 
subplot(121);
imagesc(S), axis square;
set(gca, 'clim', clim);
xlabel('Channel'), ylabel('Channel');
title('S covariance matrix');

subplot(122);
imagesc(R), axis square;
set(gca, 'clim', clim);
xlabel('Channel'), ylabel('Channel');
title('R covariance matrix');

% scree plot
% GED
[V, D] = eig(S, R);
[d, sidx] = sort(diag(D), 'descend');
V = V(:, sidx);

figure(2), clf;
plot(d, 'ks-', 'markerfacecolor', 'w', 'linew', 2, 'markersize', 13);
xlabel('Component #'), ylabel('Power ratio (\lambda)');
title('Scree plot');

% reconstruct top 2 components - should see alpha source in the static
% spectrum
cmpdat = V(:, 1:2)'*reshape(EEG.data(1:EEG.nbchan, :, :), EEG.nbchan, []);
EEG.data(EEG.nbchan+1:EEG.nbchan+2, :, :) = reshape(cmpdat, [2 EEG.pnts EEG.trials]);

comp1 = V(:, 1)'*S;
comp2 = V(:, 2)'*S;

% plot the two components
plot_simEEG(EEG, 65, 3);
subplot(221);
topoplotIndie(V(:, 1)'*S, EEG.chanlocs);
title('Component 1');
subplot(223);

plot_simEEG(EEG, 66, 4);
subplot(221);
topoplotIndie(V(:, 2)'*S, EEG.chanlocs);
title('Component 1');
subplot(223);

% filter the data or wavelet conv - extract phase and compute phase
% synchronization over time between the 2 components
% remember the data is on channel 65 and 66
frex = linspace(2, 50, 100);
waves = 2*(linspace(3, 10, length(frex))./(2*pi*frex)).^2;

wavetime = -2:1/EEG.srate:2;
halfw = floor(length(wavetime)/2)+1;
nConv = EEG.pnts * EEG.trials + length(wavetime) - 1;

dataX1 = fft(reshape(EEG.data(65, :, :), 1, []), nConv);
dataX2 = fft(reshape(EEG.data(66, :, :), 1, []), nConv);


ispc = zeros(length(frex), EEG.pnts);
pow = zeros(2, length(frex), EEG.pnts);

for fi=1:length(frex)
    waveX = fft(exp(2*1i*pi*frex(fi)*wavetime).*exp(-wavetime.^2/waves(fi)), nConv);
    waveX = waveX ./ max(waveX);

    as = ifft(waveX .* dataX1);
    as = reshape(as(halfw:end-halfw+1), [EEG.pnts, EEG.trials]);

    as2 = ifft(waveX .* dataX2);
    as2 = reshape(as2(halfw:end-halfw+1), [EEG.pnts, EEG.trials]);

    cdd = exp(1i*(angle(as) - angle(as2)));

    ispc(fi, :) = abs(mean(cdd, 2));
    pow(1, fi, :) = mean(abs(as).^2, 2);
    pow(2, fi, :) = mean(abs(as2).^2, 2);
end

figure(5), clf;
subplot(211);
plot(EEG.times, angle(as(1:EEG.pnts)), 'k', 'LineWidth', 1);
axis xy;

subplot(212);
plot(frex, mean(ispc, 2), 'ks-', 'LineWidth', 1, 'MarkerSize', 6);
axis tight;
axis xy;

figure(6), clf;
subplot(211);
plot(EEG.times, abs(as(1:EEG.pnts)).^2);
hold on;
plot(EEG.times, abs(as2(1:EEG.pnts)).^2);
xlabel('Time (ms)'), ylabel('Power (\muV^2)')
title('Example amplitude time series')
axis xy;

subplot(212);
plot(frex, mean(squeeze(pow(1, :, :)), 2), 'ks-', 'LineWidth', 1, 'color', 'blue');
hold on;
plot(frex, mean(squeeze(pow(2, :, :)), 2), 'ks-', 'LineWidth', 1, 'color', 'red');
xlabel('Freq'), ylabel('Power (\muV^2)')
title('Power spectrum')

% go back to the data and get the first 8 components (the ones without
% noise)
components8 = reshape(V(:, 1:8)' * reshape(squeeze(EEG.data(1:EEG.nbchan, :, :)), EEG.nbchan, []), 8, EEG.pnts, EEG.trials);

% compute amplitude and get the power spectrum (sum of all amplitude time
% series)


dataXs = zeros(8, nConv);

for i=1:8
    dataXs(i, :) = fft(reshape(components8(i, :, :), 1, []), nConv);
end

pow8s = zeros(8, length(frex), EEG.pnts, EEG.trials);


for fi=1:length(frex)
    waveX = fft(exp(2*1i*pi*frex(fi)*wavetime).*exp(-wavetime.^2/waves(fi)), nConv);
    waveX = waveX ./ max(waveX);

    for i=1:8
        as = ifft(waveX .* dataXs(i, :));
        as = reshape(as(halfw:end-halfw+1), [EEG.pnts, EEG.trials]);
        pow8s(i, fi, :, :) = abs(as).^2;
    end
end

% correlate amplitude (power) time series for all 8 components - compute a
% correlation matrix 8x8


% 2 figures - power time series of 8
% correlation matrix between all 8

figure(7), clf;
subplot(211);
for i=1:8
    allPow = mean(squeeze(pow8s(i, :, :)), 1);
    plot(EEG.times, allPow(1:EEG.pnts));
    hold on;
end
xlabel('Time'), ylabel('Power (\muV^2)');
title('Mean Power across all 8 components');

cormat = corrcoef(reshape(pow8s, 8, [])');

subplot(212);
imagesc(cormat), axis square;
xlabel('Component'), ylabel('Component');
title('Correlation matrix');
set(gca,'clim',[-1 1]*.7);
colorbar;

% plot all topographies for all 8 components
figure(8), clf;
for i=1:8
    subplot(2,4,i);
    topoplotIndie(V(:,i)'*S,EEG.chanlocs,'numcontour',0,'plotrad',.65);
    title([ 'Component ' num2str(i) ]);
end
