% Implementation of permutation based statistics for an EEG dataset
% Cluster correction example
% Pixel correction example
load ../data/v1_laminar.mat

npnts = size(csd, 2);
ntrials = size(csd, 3);
chan2use = [ 6 7 ];

% downsample params
times2save = -.1:.01:1.14;
times2saveidx = dsearchn(timevec', times2save');

% frequency range
min_freq = 10;
max_freq = 90;
num_frex = 35;
frex = linspace(min_freq, max_freq, num_frex);

% params for Morlet wavelets
wavtime = -2:1/srate:2-1/srate;
half_wav = (length(wavtime)-1)/2;
cycRange = [ 4 10 ];
nCycles = linspace(cycRange(1), cycRange(2), num_frex);


% FFT params
nWave = length(wavtime);
nData = npnts * ntrials;
nConv = nWave + nData - 1;

% create wavelets
cmwX = zeros(num_frex, nConv);
for fi=1:num_frex
    s = nCycles(fi) / (2*pi*frex(fi));
    cmw = exp(1i*2*pi*frex(fi).*wavtime) .* exp((-wavtime.^2) ./ (2*s^2));
    tempX = fft(cmw, nConv);
    cmwX(fi, :) = tempX ./ max(tempX);
end

dataX1 = fft(reshape(csd(chan2use(1), :, :), 1, []), nConv);
dataX2 = fft(reshape(csd(chan2use(2), :, :), 1, []), nConv);

tf = zeros(2, num_frex, length(times2save), ntrials);

for fi=1:num_frex
    % run convolution
    as1 = ifft(cmwX(fi, :) .* dataX1);
    as1 = as1(half_wav+1:end-half_wav);
    as1 = reshape(as1, npnts, ntrials);

    tf(1, fi, :, :) = abs(as1(times2saveidx, :)).^2;

    as2 = ifft(cmwX(fi, :) .* dataX2);
    as2 = as2(half_wav+1:end-half_wav);
    as2 = reshape(as2, npnts, ntrials);

    tf(2, fi, :, :) = abs(as2(times2saveidx, :)).^2;
end

% power difference for comparison
diffmap = squeeze(mean(tf(2, :, :, :), 4)) - squeeze(mean(tf(1, :, :, :), 4));

%% Visualize raw power data
clim = [0 20000];
xlim = [-.1 1];

figure(1), clf;
subplot(221);
imagesc(times2save, frex, squeeze(mean(tf(1, :, :, :), 4)));
set(gca, 'clim', clim, 'ydir', 'n', 'xlim', xlim);
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title([ 'Channel ' num2str(chan2use(1)) ]);

subplot(222);
imagesc(times2save, frex, squeeze(mean(tf(2, :, :, :), 4)));
set(gca, 'clim', clim, 'ydir', 'n', 'xlim', xlim);
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title(['Channel' num2str(chan2use(2))]);

subplot(223);
imagesc(times2save, frex, diffmap);
set(gca, 'clim', [-mean(clim) mean(clim)], 'ydir', 'n', 'xlim', xlim);
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title(['Difference: channels ' num2str(chan2use(2)) ' - ' num2str(chan2use(1))]);

%% Statistics via permutation testing
pval = 0.05;

% convert p-value to z-value
zval = abs(norminv(pval));

% nr of permutation
n_permutes = 1000;

% NULL hypothesis maps
permmaps = zeros(n_permutes, num_frex, length(times2save));

% concat maps for shuffling - 1:ntrials are from chan1 and ntrials+1:end are from channel 2
tf3d = cat(3, squeeze(tf(1, :, :, :)), squeeze(tf(2, :, :, :)));

% null hypothesis generation
for permi=1:n_permutes
    randorder = randperm(size(tf3d, 3));
    temp_tf3d = tf3d(:, :, randorder);

    % compute difference map
    permmaps(permi, :, :) = squeeze(mean(temp_tf3d(:, :, 1:ntrials), 3) - mean(temp_tf3d(:, :, ntrials+1:end), 3));
end

%% compute z - and p - values based on normalized distance to H0 distributions (per pixel)
mean_h0 = squeeze(mean(permmaps));
std_h0 = squeeze(std(permmaps));

% threshold the data
zmap = (diffmap - mean_h0) ./ std_h0;

% threshold image at p-value
zmap(abs(zmap)<zval) = 0;

%% Plot the result of the permutation testing
figure(2), clf;

subplot(221);
imagesc(times2save,frex,diffmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','nor');
title('TF map of real power values');

subplot(222);
imagesc(times2save,frex,zmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
set(gca,'clim',[-10 10],'xlim',xlim,'ydir','no');
title('Thresholded TF map of Z-values');

subplot(223);
imagesc(times2save,frex,diffmap);
hold on;
contour(times2save,frex,logical(zmap),1,'linecolor','k');
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm');
title('Power values and outlined significance regions');

%% Improve the accuracy of the prediction by implementing Cluster correction
max_cluster_sizes = zeros(1, n_permutes);

for permi=1:n_permutes
    % take each permutation map and transform to z
    threshimg = squeeze(permmaps(permi, :, :));
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold the image at p-value
    threshimg(abs(threshimg)<zval) = 0;

    % find clusters
    islands = bwconncomp(threshimg);
    if numel(islands.PixelIdxList)>0
        % count sizes
        tempclustsizes = cellfun(@length, islands.PixelIdxList);

        % store size of biggest cluster
        max_cluster_sizes(permi) = max(tempclustsizes);
    end
end

% visualize histograph of maximum cluster sizes
figure(3), clf;
hist(max_cluster_sizes, 20);
xlabel('Maximum cluster sizes'), ylabel('Number of observations');
title('Expected cluster sizes under the null hypothesis');

% find cluster threshold based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_sizes, 100-(100*pval));

% find clusters in the real thresholded zmap - if they are too small make them 0
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    if numel(islands.PixelIdxList{i})<cluster_thresh
        zmap(islands.PixelIdxList{i}) = 0;
    end
end

% plot thresholded results
figure(4), clf;
subplot(221);
imagesc(times2save,frex,diffmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title('TF power, no thresholding');
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm');

subplot(222);
imagesc(times2save,frex,diffmap);
hold on;
contour(times2save,frex,logical(zmap),1,'linecolor','k');
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title('TF power with contour');
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm');

subplot(223);
imagesc(times2save,frex,zmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title('z-map, thresholded');
set(gca,'clim',[-13 13],'xlim',xlim,'ydir','normal');

%% Another modality - Extreme Pixel Correction

% initialize matrices for cluster-based correction
max_vals = zeros(n_permutes, 2);

for permi=1:n_permutes
    % get extreme values (smallest and largest)
    temp = sort(reshape(permmaps(permi, :, :), 1, []));
    max_vals(permi, :) = temp([1 end]);
end

% visualize distributions of extreme values
figure(5), clf;
histogram(max_vals(:), 90);
xlabel('Extreme H_0 mean difference values');
ylabel('Count');

% threshold the data
pthresh = .045;

tmp = sort(max_vals(:, 1));
thresh_lo = tmp(round(pthresh*n_permutes));

% repeat for upper bound
tmp = sort(max_vals(:, 2));
thresh_hi = tmp(round((1-pthresh)*n_permutes));

% plot on histogram
figure(5), hold on;
plot([1 1]*thresh_lo, get(gca, 'ylim'), 'r--', 'linew', 2);
plot([1 1]*thresh_hi, get(gca, 'ylim'), 'r--', 'linew', 2);

% threshold real data
zmap = diffmap;
zmap(zmap > thresh_lo & zmap < thresh_hi) = 0;

figure(6), clf;
subplot(221);
imagesc(times2save,frex,diffmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title('TF power map, no thresholding');
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','n');

subplot(222), hold on;
imagesc(times2save,frex,diffmap);
contour(times2save,frex,logical(zmap),1,'linecolor','k','linewidth',3);
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title('TF power map with contour');
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','normal');

subplot(223);
imagesc(times2save,frex,zmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title('TF power map, thresholded');
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','no');
