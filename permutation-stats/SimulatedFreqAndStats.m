% Example of simulated EEG data - with frequency analysis and permutation statistics

% generic params
npnts = 2501;
time = linspace(-1, 3, npnts);
srate = 1./mean(diff(time));
ntrials = 100;
nPerms = 1000;
pval = .05;

% convert p-value to z-value (note: use 1.96 for p<.05 2-tailed without stats-toolbox)
sigThresh = norminv(1-pval/2);

% some parameters copied from the solution with how the data should look
% like
freq1a = 6;
freq1b = 15;
freq2a = 7;
freq2b = 14;
peak1a = 1;
peak1b = 1.6;
peak2a = .9;
peak2b = 1.7;
amp1a  = 4;
amp1b  = 3;
amp2a  = 2;
amp2b  = 5;
fwhm1a = .6;
fwhm1b = .4;
fwhm2a = .6;
fwhm2b = .5;

% gausians
gausw1a = exp(-4 * log(2)*(time-peak1a).^2 / fwhm1a.^2);
gausw1b = exp(-4 * log(2)*(time-peak1b).^2 / fwhm1b.^2);

gausw2a = exp(-4 * log(2)*(time-peak2a).^2 / fwhm2a.^2);
gausw2b = exp(-4 * log(2)*(time-peak2b).^2 / fwhm2b.^2);

[data1, data2] = deal(zeros(npnts, ntrials));
noisestd = 1.5;

%%% noise:
% 1/f, trial-unique
ed = 50; % exponential decay parameter

for triali=1:ntrials
    sig1a = amp1a * sin(2*pi*freq1a*time + 2*pi*randn) .* gausw1a;
    sig1b = amp1b * sin(2*pi*freq1b*time + 2*pi*randn) .* gausw1b;

    % generate 1/f noise
    as = rand(1,floor(npnts/2)-1) .* exp(-(1:floor(npnts/2)-1)/ed);
    as = [as(1) as 0 0 as(:,end:-1:1)];
    fc = as .* exp(1i*2*pi*rand(size(as)));
    noise = real(ifft(fc));
    noise = noise*(noisestd/std(noise));

    data1(:, triali) = sig1a + sig1b + noise;

    sig2a = amp2a * sin(2*pi*freq2a*time + 2*pi*randn) .* gausw2a;
    sig2b = amp2b * sin(2*pi*freq2b*time + 2*pi*randn) .* gausw2b;

    as = rand(1,floor(npnts/2)-1) .* exp(-(1:floor(npnts/2)-1)/ed);
    as = [as(1) as 0 0 as(:,end:-1:1)];
    fc = as .* exp(1i*2*pi*rand(size(as)));
    noise = real(ifft(fc));
    noise = noise*(noisestd/std(noise));

    data2(:, triali) = sig2a + sig2b + noise;
end

% plot the two signals
figure(1), clf;
subplot(211);
plot(time, data1(:, 1), time, data2(:, 1));
ylabel('Amplitude');
xlabel('Time (ms)');
legend({'data1 - trial 1', 'data2 - trial2'});
title('Signal and noise example on trial 1');

subplot(212);
plot(time, mean(data1, 2), time, mean(data2, 2));
ylabel('Amplitude');
xlabel('Time (ms)');
legend({'data1', 'data2'});
title('Signal and noise mean all trials');

%%% data analysis:
% wavelet convolution

% parameters for wavelet convolution
num_frex = 60;
min_frex = 2;
max_frex = 30;
frex = linspace(min_frex, max_frex, num_frex);
wavtime = -2:1/srate:2;
half_wav = (length(wavtime)-1)/2;

nWave = length(wavtime);
nData = npnts * ntrials;
nConv = nWave + nData - 1;

cmwX = zeros(num_frex, nConv);
fwhm = linspace(.6,1,num_frex);

% generate Morlet Wavelets
for fi=1:num_frex
    % create wavelet
    cmw = exp(1i*2*pi*frex(fi).*wavtime).*exp(-4*log(2)*wavtime.^2 / fwhm(fi)^2);
    tempX = fft(cmw, nConv);
    cmwX(fi, :) = tempX ./ max(tempX);
end

tf = zeros(2, num_frex, npnts, ntrials);
dataX1 = fft(reshape(data1, 1, []), nConv);
dataX2 = fft(reshape(data2, 1, []), nConv);

for fi=1:num_frex
    % convolution
    as1 = ifft(cmwX(fi, :) .* dataX1);
    as1 = as1(half_wav+1:end-half_wav);
    as1 = reshape(as1, npnts, ntrials);
    tf(1, fi, :, :) = abs(as1).^2;

    as2 = ifft(cmwX(fi, :) .* dataX2);
    as2 = as2(half_wav+1:end-half_wav);
    as2 = reshape(as2, npnts, ntrials);
    tf(2, fi, :, :) = abs(as2).^2;
end

figure(2), clf;
subplot(211)
imagesc(time, frex, mean(squeeze(tf(1, :, :, :)), 3));
axis xy
set(gca,'clim',[0 .22])
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Condition1');

subplot(212)
imagesc(time, frex, mean(squeeze(tf(2, :, :, :)), 3));
axis xy
set(gca,'clim',[0 .22])
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Condition1');


%%% statistics
% power difference
diffmap = squeeze(mean(tf(2, :, :, :), 4)) - squeeze(mean(tf(1, :, :, :), 4));

% visualize diff map
figure(3), clf;
subplot(221);
contourf(time, frex, diffmap, 40, 'linecolor', 'none');
set(gca, 'clim', [-.9 1], 'ydir', 'n', 'xlim', [-1 3]);
xlabel('Time (ms)'), ylabel('Frequenccy (Hz)');
title(['Diff map between condition 2 and 1']);

% mapwise permutation testing
pval = 0.05;
zval = abs(norminv(pval));
n_permutes = 100;
permmaps = zeros(n_permutes, num_frex, npnts);

% concat maps for shuffling
tf3d = cat(3, squeeze(tf(1, :, :, :)), squeeze(tf(2, :, :, :)));
condlabels = (1:ntrials*2)>ntrials;

% null hypothesis generation
for permi=1:n_permutes
    randorder = randperm(size(tf3d, 3));
    temp_tf3d = tf3d(:, :, randorder);

    permmaps(permi, :, :) = squeeze(mean(temp_tf3d(:, :, 1:ntrials), 3) - mean(temp_tf3d(:, :, ntrials+1:end), 3));
end

% compute z and p values based on normalized distance to H0 distributions
mean_h0 = squeeze(mean(permmaps));
std_h0 = squeeze(std(permmaps));

% threshold the data
zmap = (diffmap - mean_h0) ./ std_h0;

thresh_zmap = zmap;
thresh_zmap(abs(thresh_zmap)<zval) = 0;

% visualization
subplot(222);
contourf(time, frex, zmap, 40, 'linecolor', 'none');
set(gca, 'clim', [-7 7], 'ydir', 'n', 'xlim', [-1 3]);
xlabel('Time (ms)'), ylabel('Frequenccy (Hz)');
title('Difference z-map');

subplot(223);
contourf(time, frex, thresh_zmap, 40, 'linecolor', 'none');
set(gca, 'clim', [-7 7], 'ydir', 'n', 'xlim', [-1 3]);
xlabel('Time (ms)'), ylabel('Frequenccy (Hz)');
title('Thresholded difference z-map');

subplot(224);
contourf(time, frex, diffmap, 40, 'linecolor', 'none');
hold on;
contour(time,frex,logical(thresh_zmap),1,'linecolor','k');
set(gca, 'clim', [-.9 1], 'ydir', 'n', 'xlim', [-1 3]);
xlabel('Time (ms)'), ylabel('Frequenccy (Hz)');
title('Difference map with thresholded outlines');

% cluster correction
% separately for positive vs. negative features!
positive_max_cluster_sizes = zeros(1, n_permutes);
negative_max_cluster_sizes = zeros(1, n_permutes);

for permi=1:n_permutes
    threshimg = squeeze(permmaps(permi, :, :));
    threshimg = (threshimg - mean_h0) ./ std_h0;

    %threshold the image at p-value
    positive_thresh = threshimg;
    positive_thresh(positive_thresh > 0 & abs(positive_thresh - mean_h0)<zval) = 0;
    positive_thresh(positive_thresh < 0) = 0;

    negative_thresh = threshimg;
    negative_thresh(negative_thresh < 0 & abs(negative_thresh - mean_h0)<zval) = 0;
    negative_thresh(negative_thresh > 0) = 0;

    islands = bwconncomp(positive_thresh);
    if numel(islands.PixelIdxList)>0
        tempclustsizes = cellfun(@length, islands.PixelIdxList);
        positive_max_cluster_sizes(permi) = max(tempclustsizes);
    end

    islands = bwconncomp(negative_thresh);
    if numel(islands.PixelIdxList)>0
        tempclustsizes = cellfun(@length, islands.PixelIdxList);
        negative_max_cluster_sizes(permi) = max(tempclustsizes);
    end
end

figure(4), clf;
subplot(211);
hist(positive_max_cluster_sizes, 20);
xlabel('Maximum cluster sizes'), ylabel('Number of observations');
title('Positive Clusters');

subplot(212);
hist(negative_max_cluster_sizes, 20);
xlabel('Maximum cluster sizes'), ylabel('Number of observations');
title('Positive Clusters');

positive_cluster_thresh = prctile(positive_max_cluster_sizes, 100-(100*pval));
negative_cluster_thresh = prctile(negative_max_cluster_sizes, 100-(100*pval));

cluster_correct_zmap = thresh_zmap;
islands = bwconncomp(cluster_correct_zmap);
for i=1:islands.NumObjects
    island_size = numel(islands.PixelIdxList{i});
    island_mean_value = mean(zmap(islands.PixelIdxList{i}));

    % check for positive
    if island_mean_value > 0 && (island_size < positive_cluster_thresh)
        cluster_correct_zmap(islands.PixelIdxList{i}) = 0;
    end

    % check for negative
    if island_mean_value < 0 && (island_size < positive_cluster_thresh)
        cluster_correct_zmap(islands.PixelIdxList{i}) = 0;
    end
end

figure(5), clf;
contourf(time, frex, thresh_zmap, 40, 'linecolor', 'none');
hold on;
contour(time,frex,logical(cluster_correct_zmap),1,'linecolor','k');
set(gca, 'clim', [-7 7], 'ydir', 'n', 'xlim', [-1 3]);
xlabel('Time (ms)'), ylabel('Frequenccy (Hz)');
title('Difference map with thresholded outlines, cluster corrected');
