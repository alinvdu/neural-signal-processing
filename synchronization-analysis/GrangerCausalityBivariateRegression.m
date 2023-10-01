% Directed synchrony through Granger Causality (Prediction)

load ../data/sampleEEGdata.mat

% define channels to compute granger synchrony between
chan1name = 'fcz';
chan2name = 'o1';

% find de indices of these channels
chan1 = find(strcmpi(chan1name, {EEG.chanlocs.labels}));
chan2 = find(strcmpi(chan2name, {EEG.chanlocs.labels}));

% define autoregression params
order = 14;

% get autoregression coeefcs andd error for each signal
[Ax, Ex] = armorf(EEG.data(chan1, :, 1), 1, EEG.pnts, order);
[Ay, Ey] = armorf(EEG.data(chan2, :, 1), 1, EEG.pnts, order);

% Bivariate autogression and associated error term
[Axy, E] = armorf(EEG.data([chan1 chan2], :, 1), 1, EEG.pnts, order);

% time-domain causal estimate
granger_chan2_to_chan1 = log(Ex/E(1, 1));
granger_chan1_to_chan2 = log(Ey/E(2, 2));

disp(['Granger prediction from ' chan1name ' to ' chan2name ' is ' num2str(granger_chan1_to_chan2)]);
disp(['Granger prediction from ' chan2name ' to ' chan1name ' is ' num2str(granger_chan2_to_chan1)]);

%% compute over time windows
x2yT = zeros(1, EEG.pnts);
y2xT = zeros(1, EEG.pnts);

% GC params
iwin = 300;
iorder = 15;

% convert window/order to points
win = round(iwin/(1000/EEG.srate));
order = round(iorder/(1000/EEG.srate));

for timei=1:EEG.pnts-win
    % data from all trials in this time window andd normalized
    tempdata = zscore(reshape(EEG.data([chan1 chan2], timei:timei+win-1, 1), 2, win), 0, 2);

    % fit AR moddel
    [Ax, Ex] = armorf(tempdata(1, :), 1, win, order);
    [Ay, Ey] = armorf(tempdata(2, :), 1, win, order);
    [Axy, E] = armorf(tempdata, 1, win, order);

    % causal estimate
    y2xT(timei) = log(Ex/E(1, 1));
    x2yT(timei) = log(Ey/E(2, 2));
end

figure(1), clf, hold on
plot(EEG.times, x2yT);
plot(EEG.times, y2xT, 'r');
legend({[ 'GC: ' chan1name ' -> ' chan2name ]; [ 'GC' chan2name ' -> ' chan1name ]});

title(['Window length: ' num2str(win) ' ms, order: ' num2str(iorder) ' ms' ]);
xlabel('Time (ms)');
ylabel('Granger prediction estimate');
set(gca, 'xlim', [-200 200]);
