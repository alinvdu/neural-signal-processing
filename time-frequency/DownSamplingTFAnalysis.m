clear
load v1_laminar

% the time points we want to save (the downsampling times vector - way smaller than times)
times2save = -.2:.025:1.2;

% find indices
tidx = dsearchn(timevec', times2save');

freqrange = [10 100];
numfrex = 42;

% set up convolution parameters
wavtime = -1:1/srate:1-1/srate;
frex = linspace(freqrange(1), freqrange(2), numfrex);
nData = size(csd, 2) * size(csd, 3);
nKern = length(wavtime);
nConv = nData + nKern - 1;
halfwav = (length(wavtime)-1)/2;

% create wavelets
cmwX = zeros(numfrex, nConv);
for fi=1:numfrex
    % create time-domain complex Morlet Wavelet
    cmw = exp(1i*2*pi*frex(fi).*wavtime).*exp(-4*log(2)*wavtime.^2 / .3^2);

    % compute fourier coefficients of wavelet and normalize
    cmwX(fi, :) = fft(cmw, nConv);
    cmwX(fi, :) = cmwX(fi, :) ./ max(cmwX(fi, :));
end

% initialize time-frequency output matrix
tf = cell(2, 1);

% compute Fourier coefficients of EEG data
eegX = fft(reshape(csd(7, :, :), 1, []), nConv);

for fi=1:numfrex
    % second and third steps of convolution
    as = ifft(eegX .* cmwX(fi, :));

    % cut wavelet back to size of data
    as = as(halfwav+1:end-halfwav);
    as = reshape(as, length(timevec), size(csd, 3));

    tf{1}(fi, :) = mean(abs(as).^2, 2);
    tf{2}(fi, :) = mean(abs(as(tidx, :)).^2, 2);
end

% visualized downsampled data
figure(18), clf
titles = {'Full'; 'Downsampled'};

for i=1:2
    subplot(1, 2, i);

    if i==1, t=timevec; else, t=times2save; end
    
    contourf(t, frex, tf{i}, 40, 'linecolor', 'none');
    set(gca, 'clim', [0 10000], 'xlim', [-.2 1.2]);
    xlabel('Time (s)'), ylabel('Frequencies (Hz)');
    title([titles{i} 'Time-Frequency power'])
end
