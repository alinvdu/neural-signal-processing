load emptyEEG

% dipole location
diploc1 = 109;
diploc2 = 118;

% plot brain dipoles
figure(1), clf, subplot(131)
plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), 'bo','markerfacecolor','y')
hold on
plot3(lf.GridLoc(diploc1,1), lf.GridLoc(diploc1,2), lf.GridLoc(diploc1,3), 'ks','markerfacecolor','k','markersize',10)
plot3(lf.GridLoc(diploc2,1), lf.GridLoc(diploc2,2), lf.GridLoc(diploc2,3), 'rs','markerfacecolor','r','markersize',10)
rotate3d on, axis square
title('Brain dipole locations')

% Project dipole
subplot(132)
topoplotIndie(-lf.Gain(:, 1, diploc1), EEG.chanlocs, 'numcontour', 0, 'electrodes', 'numbers','shading','interp');
set(gca,'clim',[-1 1]*40)
title('Signal dipole projection')

subplot(133)
topoplotIndie(-lf.Gain(:,1,diploc2), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[-1 1]*40)
title('Signal dipole projection')

%% EEG Params
EEG.pnts = 1143;
EEG.trials = 150;
EEG.times = (0:EEG.pnts-1)/EEG.srate - .2;

% initialize EEG data
EEG.data = zeros([ EEG.nbchan EEG.pnts EEG.trials ]);

%% create simulated data
% Gaussian
peaktime = .5; %seconds
fwhm = .12;

% create a Gaussian taper
gaus = exp( -(4*log(2)*(EEG.times-peaktime).^2) / fwhm^2 );

sineFreq1 = 9;
sineFreq2 = 14;

sine1 = sin(2 * pi * sineFreq1 * EEG.times);
sine2 = sin(2 * pi * sineFreq2 * EEG.times);

for triali=1:EEG.trials
    
    % initialize all dipole data
    dipdat = .01 * randn(size(lf.Gain,3),EEG.pnts);
    
    dipdat(diploc1,:) = sine1 .* gaus;
    dipdat(diploc2,:) = sine2 .* gaus;
    
    % compute one trial
    EEG.data(:,:,triali) = squeeze(lf.Gain(:,1,:))*dipdat;
end

%% Separate two sources based on power spectra.
% FFT over all channels
dataX = fft(EEG.data, [], 2);
dataPow = mean(abs(dataX), 3);

% freq
hz = linspace(0, EEG.srate/2, floor(EEG.pnts/2)+1);

% frequency cutoffs
frex = [4 11 20];
fidx = dsearchn(hz', frex');

% power in the first spectral window
powLow = mean(dataPow(:, fidx(1):fidx(2)), 2);
powHi = mean(dataPow(:, fidx(2):fidx(3)), 2);

figure(5), clf
subplot(221)
topoplotIndie(-lf.Gain(:,1,diploc1), EEG.chanlocs,'numcontour',0,'electrodes','off','shading','interp');
title('Ground truth, dipole 1')

subplot(222)
topoplotIndie(-lf.Gain(:,1,diploc2), EEG.chanlocs,'numcontour',0,'electrodes','off','shading','interp');
title('Ground truth, dipole 2')

subplot(223)
topoplotIndie(powLo,EEG.chanlocs,'numcontour',0,'electrodes','off','shading','interp');
title([ 'Data: ' num2str(round(hz(fidx(1)),2)) ' - ' num2str(round(hz(fidx(2)),2)) ' Hz' ])

subplot(224)
topoplotIndie(powHi,EEG.chanlocs,'numcontour',0,'electrodes','off','shading','interp');
title([ 'Data: ' num2str(round(hz(fidx(2)),2)) ' - ' num2str(round(hz(fidx(3)),2)) ' Hz' ])
