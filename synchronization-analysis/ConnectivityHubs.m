% connectivity hubs - regions of powerful connections in the brain
load ../data/sampleEEGdata.mat
EEG.data = double(EEG.data); % it's a good idea to do this for analysis

frex = 10;
tidx = dsearchn(EEG.times', [0 500]');

% wavelet params
wtime = -1:1/EEG.srate:1;
fwhm = .1;

% convolution params
nData = EEG.pnts * EEG.trials;
nWave = length(wtime);
nConv = nData + nWave - 1;
halfW = floor(nWave/2);

% create wavelet
cmw = fft(exp(1i*2*pi*frex*wtime).*exp(-4*log(2)*wtime.^2 / fwhm^2), nConv);

% fourier spectrum for data
dataX = fft(reshape(EEG.data, EEG.nbchan, []), nConv, 2);

% convolution. All channels at once
as = ifft(dataX .* cmw);
as = as(:, halfW:end-halfW-1);
as = reshape(as, size(EEG.data));

allphases = angle(as);

% compute all to all PLI
pliall = zeros(EEG.nbchan);

for chani=1:EEG.nbchan
    for chanj=chani+1:EEG.nbchan
        % euler-format phase diffs
        cdd = exp(1i*(allphases(chani, tidx(1):tidx(2), :)-allphases(chanj, tidx(1):tidx(2), :)));

        % compute PLI for this channel pair
        plitmp = mean(abs(mean(sign(imag(cdd)), 2)));

        pliall(chani, chanj) = plitmp;
        pliall(chanj, chani) = plitmp; % symmetric
    end
end

% visualize matrix of connections
figure(1), clf;
imagesc(pliall);
axis square;
xlabel('Channels'), ylabel('Channels');
title(['All-to-all connectivity at ' num2str(frex) ' Hz']);
set(gca, 'clim', [.2 .6]);

% calculate hubness
% define thresholdd
distdata = nonzeros(triu(pliall));

% define a threshold
thresh = median(distdata) + std(distdata);

% binarize the matrix by testing against the threshold
pliallThresh = pliall > thresh;

% plot results
figure(2), clf;
subplot(311), hold on;
histogram(distdata, 50); % for the threshold
plot([1 1]*thresh, get(gca, 'ylim'), 'r--', 'linew', 3)
xlabel('PLI (synch. strength)'), ylabel('Count');
legend({'Distribution'; 'Threshold'});

subplot(3, 2, [3 5]);
imagesc(pliallThresh);
axis square;
xlabel('Channels'), ylabel('Channels');
title([ 'All-to-all connectivity at ' num2str(frex) ' Hz' ]);

subplot(3, 2, [4 6]);
topoplotIndie(sum(pliallThresh)/(EEG.nbchan-1), EEG.chanlocs, 'numcontour', 0);
set(gca, 'clim', [.1 .4]);
title('Topoplot of "hubness"');
colormap hot;
colorbar;
