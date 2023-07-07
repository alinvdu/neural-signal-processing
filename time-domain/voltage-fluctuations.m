load ../data/sampleEEGdata.mat

% channel to pick
chan2use = 'o1';

% time window for negative peak
% Mike's choices: 
negpeaktime = [  50 110 ];
pospeaktime = [ 110 170 ];
indicesneg = arrayfun(@(x) find(abs(EEG.times(:) - x) == min(abs(EEG.times(:) - x)), 1), negpeaktime);
indicespos = arrayfun(@(x) find(abs(EEG.times(:) - x) == min(abs(EEG.times(:) - x)), 1), pospeaktime);

% compute channel index
chanidx = find(strcmpi({EEG.chanlocs.labels}, chan2use));

%%% compute ERP
erp = double( mean(EEG.data(chanidx,:,:),3) );

% plot ERP
figure(1), clf
plot(EEG.times,erp,'k','linew',1)
set(gca,'xlim',[-300 1000])

% plot patches over areas
ylim = get(gca,'ylim');
ph = patch(negpeaktime([1 1 2 2]),ylim([1 2 2 1]),'y');
set(ph,'facealpha',.8,'facecolor','yellow', 'edgecolor', 'none')

ph2 = patch(pospeaktime([1 1 2 2]),ylim([1 2 2 1]),'y');
set(ph2,'facealpha',.8,'facecolor','green', 'edgecolor', 'none')

% move the patches to the background
set(gca,'Children',flipud( get(gca,'Children') ))

%% first low-pass filter (windowed sinc function)

lowcut = 15;
filttime = -.3:1/EEG.srate:.3;
filtkern = sin(2*pi*lowcut*filttime) ./ filttime;

% adjust NaN and normalize filter to unit-gain
filtkern(~isfinite(filtkern)) = max(filtkern);
filtkern = filtkern./sum(filtkern);

% windowed sinc filter
filtkern = filtkern .* hann(length(filttime))';

% inspect the filter kernel
figure(2), clf
subplot(211)
plot(filttime, filtkern,'k','linew',2)
xlabel('Time (s)')
title('Time domain')


subplot(212)
hz = linspace(0,EEG.srate,length(filtkern));
plot(hz, abs(fft(filtkern)).^2,'ks-','linew',2)
set(gca,'xlim',[0 lowcut*3])
xlabel('Frequency (Hz)'), ylabel('Gain')
title('Frequency domain')

%% now filter the ERP and replot

% apply filter
filteredErp = filtfilt(filtkern, 1, erp);

% plot on top of unfiltered ERP
figure(1)
hold on
plot(EEG.times,filteredErp,'k','linew',2, 'color', 'red')

%% peak-to-peak voltages and timings

%%%% first for unfiltered ERP

% find minimum/maximum peak values and peak times
minimumerp = erp(indicesneg(1):indicesneg(2));
minimum = min(minimumerp);

maximumerp = erp(indicespos(1):indicespos(2));
maximum = max(maximumerp);

minimumIndexErp = find(erp == minimum);
maximumIndexErp = find(erp == maximum);

% get results (peak-to-peak voltage and latency)
erpP2P = abs(maximum - minimum);
erpP2Plat = EEG.times(maximumIndexErp) - EEG.times(minimumIndexErp);


%%%% then for low-pass filtered ERP

% find minimum/maximum peak values and peak times
minimumerp = filteredErp(indicesneg(1):indicesneg(2));
minimum = min(minimumerp);

maximumerp = filteredErp(indicespos(1):indicespos(2));
maximum = max(maximumerp);

minimumIndexUn = find(filteredErp == minimum);
maximumIndexUn = find(filteredErp == maximum);


% get results (peak-to-peak voltage and latency)
erpFP2P = abs(maximum - minimum);
erpFP2Plat = EEG.times(maximumIndexUn) - EEG.times(minimumIndexUn);

%% Report the results in the command window

% clear the screen
clc

fprintf('\nRESULTS FOR PEAK POINT:')
fprintf('\n   Peak-to-peak on unfiltered ERP: %5.4g muV, %4.3g ms span.',erpP2P,erpP2Plat)
fprintf('\n   Peak-to-peak on filtered ERP:   %5.4g muV, %4.3g ms span.\n\n',erpFP2P,erpFP2Plat)

%%

%% repeat for mean around the peak

% time window for averaging (one-sided!!)
win = 10; % in ms
% now convert to indices
eegTimeMin = EEG.times(minimumIndexErp);
minTimeValues = [eegTimeMin - 10 eegTimeMin + 10];
minIndicesValues = arrayfun(@(x) find(abs(EEG.times(:) - x) == min(abs(EEG.times(:) - x)), 1), minTimeValues);
meanMined = mean(erp(minIndicesValues(1):minIndicesValues(2)));

eegTimeMax = EEG.times(maximumIndexErp);
maxTimeValues = [eegTimeMax - 10  eegTimeMax + 10];
maxIndicesValues = arrayfun(@(x) find(abs(EEG.times(:) - x) == min(abs(EEG.times(:) - x)), 1), maxTimeValues);
meanMaxed = mean(erp(maxIndicesValues(1):maxIndicesValues(2)));

% get results (peak-to-peak voltage and latency)
erpP2P = meanMaxed - meanMined;

% now for unfiltered
eegTimeMin = EEG.times(minimumIndexUn);
minTimeValues = [eegTimeMin - 10 eegTimeMin + 10];
minIndicesValues = arrayfun(@(x) find(abs(EEG.times(:) - x) == min(abs(EEG.times(:) - x)), 1), minTimeValues);
meanMined = mean(filteredErp(minIndicesValues(1):minIndicesValues(2)));

eegTimeMax = EEG.times(maximumIndexUn);
maxTimeValues = [eegTimeMax - 10  eegTimeMax + 10];
maxIndicesValues = arrayfun(@(x) find(abs(EEG.times(:) - x) == min(abs(EEG.times(:) - x)), 1), maxTimeValues);
meanMaxed = mean(filteredErp(maxIndicesValues(1):maxIndicesValues(2)));

% get results (peak-to-peak voltage and latency)
erpFP2P = meanMaxed - meanMined;


%% Report the results in the command window

fprintf('\nRESULTS FOR WINDOW AROUND PEAK:')
fprintf('\n   Peak-to-peak on unfiltered ERP: %5.4g muV, %4.3g ms span.',erpP2P,erpP2Plat)
fprintf('\n   Peak-to-peak on filtered ERP:   %5.4g muV, %4.3g ms span.\n\n',erpFP2P,erpFP2Plat)

%% done.
