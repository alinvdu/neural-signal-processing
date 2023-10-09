clear
load ../data/groupTFdata.mat

nSubs = size(tf, 1);

% plot the mean of the 2 conditions
figure(1), clf;
for i=1:2
    subplot(1, 2, i);
    contourf(times2save, frex, squeeze(mean(tf(:, i, :, :))), 40, 'linecolor', 'none');
    set(gca, 'clim', [-1 1]*3);
    xlabel('Time (ms)'), ylabel('Frequency (Hz)');
    title(['Condition' num2str(i) ]);
    axis square;
end

% plot each subject for condition 2
figure(2), clf;
for i=1:25
    subplot(5, 5, i);
    contourf(times2save, frex, squeeze(tf(i, 2, :, :)), 40, 'linecolor', 'none');
    set(gca, 'clim', [-1 1]*3, 'xtick', [], 'ytick', []);
    axis square;
    title([ 's' num2str(i) ]);
end

% total average
grandAve = squeeze(mean(mean(tf, 2), 1));
figure(3), clf;
contourf(times2save,frex,grandAve,40,'linecolor','none');
set(gca,'clim',[-1 1]*2,'YScale','lo','ytick',ceil(frex(1:8:end)));
axis square;
xlabel('Time (ms)'), ylabel('Frequency (Hz)');
title('Total average over subjects and conditions');

%% pick time-frequency regions based on condition average
time1 = [285 630];
time2 = [180 555];

freq1 = [3 7];
freq2 = [15 30];

% draw boxes on top
hold on
rectangle('Position',[time1(1) freq1(1) diff(time1) diff(freq1)],'linew',3);
rectangle('Position',[time2(1) freq2(1) diff(time2) diff(freq2)],'linew',3);
text(time1(1),freq1(2),'ROI 1','VerticalAlignment','bottom','FontSize',20);
text(time2(1),freq2(2),'ROI 2','VerticalAlignment','bottom','FontSize',20);

% Extract data for futher analysis
time1idx = dsearchn(times2save',time1');
time2idx = dsearchn(times2save',time2');

freq1idx = dsearchn(frex',freq1');
freq2idx = dsearchn(frex',freq2');

labels = { 'theta cond1'; 'theta cond2'; 'beta cond1'; 'beta cond2' };
data = zeros(4, nSubs);

data(1,:) = squeeze(mean(mean(tf(:, 1, freq1idx(1):freq1idx(2), time1idx(1):time1idx(2)), 3), 4));
data(2,:) = squeeze(mean(mean(tf(:, 2, freq1idx(1):freq1idx(2), time1idx(1):time1idx(2)), 3), 4));
data(3,:) = squeeze(mean(mean(tf(:, 1, freq2idx(1):freq2idx(2), time2idx(1):time2idx(2)), 3), 4));
data(2,:) = squeeze(mean(mean(tf(:, 2, freq2idx(1):freq2idx(2), time2idx(1):time2idx(2)), 3), 4));

figure(4), clf, hold on;
bh = bar([1 2 4 5],mean(data,2));
eh = errorbar([1 2 4 5],mean(data,2),std(data,[],2)/sqrt(nSubs-1),'.');
set(gca,'XTickLabel',labels,'xtick',[1 2 4 5]);
set(bh,'FaceColor',[.7 .3 .9]);
set(eh,'LineWidth',10,'color','k');
ylabel('Power (dB)');