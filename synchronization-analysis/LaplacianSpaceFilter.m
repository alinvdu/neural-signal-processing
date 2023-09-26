% Surface Laplacian Spatial Filtering - estimates the second spatial derivative
% in two dimensions across scalp features. This method is aimed at enhancing the
% spatial resolution of the signals and reducing the effects of volume conduction.

clear
load ../data/sampleEEGdata.mat

time2plot = 250;
chan2plot = 'cz';

% convert to indices
tidx = dsearchn(EEG.times', time2plot);
chanidx = strcmpi({EEG.chanlocs.labels}, chan2plot);

% compute laplacian and store as a new field
EEG.lap = laplacian_perrinX(EEG.data, [EEG.chanlocs.X], [EEG.chanlocs.Y], [EEG.chanlocs.Z]);

% The voltage and Laplacian are in different scales, z normalize.
voltERP = mean(EEG.data(chanidx, :, :), 3);
voltERP = zscore(voltERP);

lapERP = mean(EEG.lap(chanidx, :, :), 3);
lapERP = zscore(lapERP);

% plotting
figure(1), clf;
subplot(221);
topoplotIndie(mean(EEG.data(:, tidx, :), 3), EEG.chanlocs, 'electroddes', 'labels', 'numcontour', 0);
title(['Voltage ()' num2str(time2plot) ' ms)']);

subplot(222);
topoplotIndie(mean(EEG.lap(:, tidx, :), 3), EEG.chanlocs, 'electrodes', 'numbers', 'numcontour', 0);
title(['Laplacian (' num2str(time2plot) ' ms)']);

subplot(212);
plot(EEG.times, voltERP, EEG.times, lapERP, 'linew', 2);
set(gca, 'xlim', [-300 1200]);
legend({'Voltage'; 'Laplacian'});
title(['ERP from channel' chan2plot]);
xlabel('Time (ms)'), ylabel('Data (z-score)');
