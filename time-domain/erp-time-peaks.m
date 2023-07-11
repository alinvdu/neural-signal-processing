load ../data/sampleEEGdata.mat

times = [ 100 400 ];
times_ind = dsearchn(EEG.times', times')';

peak_times = zeros(length(EEG.times));

%%
for chanloci=1:length(EEG.chanlocs)
    erp = double(mean(EEG.data(chanloci, times_ind(1):times_ind(2), :), 3));
    [v, index] = max(erp);
    peak_times(chanloci) = EEG.times(times_ind(1) + index - 1);
end

figure(1)
subplot(2, 1, 1);
topoplotIndie(peak_times,EEG.chanlocs, 'numcontour', 4, 'electrodes', 'numbers');
title('Unfiltered ERP peak times (100-400ms)');
set(gca, 'clim', times);
colormap hot
colorbar

%% first low-pass filter (windowed sinc function)
lowcut = 15;
filttime = -.3:1/EEG.srate:.3;
filtkern = sin(2*pi*lowcut*filttime) ./ filttime;

% adjust NaN and normalize filter to unit-gain
filtkern(~isfinite(filtkern)) = max(filtkern);
filtkern = filtkern./sum(filtkern);

% windowed sinc filter
filtkern = filtkern .* hann(length(filttime))';

peak_times_filtered = zeros(length(EEG.times));

%%
for chanloci=1:length(EEG.chanlocs)
    erp = double(mean(EEG.data(chanloci, :, :), 3));
    filtered_erp = filtfilt(filtkern, 1, erp);
    filtered_erp_sliced = filtered_erp((times_ind(1):times_ind(2)));
    [v, index] = max(filtered_erp_sliced);
    peak_times_filtered(chanloci) = EEG.times(times_ind(1) + index - 1);
end

subplot(2, 1, 2);
topoplotIndie(peak_times_filtered,EEG.chanlocs, 'numcontour', 4, 'electrodes', 'numbers');
title('Filtered ERP peak times (100-400ms)');
set(gca, 'clim', times);
colormap hot
colorbar