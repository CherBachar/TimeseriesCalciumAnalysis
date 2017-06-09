function [handles] = detectPeaks(trace, traceNeuropil, handles)

Plot = 0;
validMask = [1 1 1 0 -1 -1 -1];
numCells = size(trace,1);
time = size(trace,2);
thresh = mean(trace(:))+std(trace(:));
sizeImage = handles.sizeImage;
ITimeseries = handles.ITimeseries;
Cells = handles.Cells;
sumPixels = handles.sumPixels;
tau1 = 10;
tau2 = 10;

%% df/d calculation
sort_trace_by_F = sort(trace, 2, 'ascend');
for x = 1:numCells
    min_10(x, 1) = mean(sort_trace_by_F(x, 1:10));
end

% (F-F0) / F0
df_fixedF0WOBack = zeros(numCells,time);
for c=1:numCells
    for t=1:(time-1)
        df_fixedF0WOBack(c,t) = (trace(c,t+1)- min_10(c))/(min_10(c));
    end
    df_fixedF0WOBack(c,t) = (trace(c,t+1)- min_10(c))/(min_10(c));
end
df_fixedF0WOBack(:,time) = df_fixedF0WOBack(:,time-1);


%smoothing
for i = 1:numCells
    df_fixedF0(i,:) = smooth(df_fixedF0WOBack(i,:), 'loess');
end


%% df/d calculation using the neuropil (i.e. surrounding area) to normalize

%with smoothing
for i = 1:numCells
    traceNeuropilSmooth(i,:) = smooth(traceNeuropil(i,:), 'loess');
end
min10mat = repmat(min_10, [1, time]);
df_fixedF0Back = ((trace - traceNeuropilSmooth)-min10mat)./min10mat;


%% df/f calculation from papaer: In vivo two-photon imaging of sensory-evoked dendritic calcium signals in cortical neurons
%1.Calculate the mean fluorescence F(t) of a region of interest (ROI) for each time point t
Ft = zeros(numCells, time);
for c = 1:numCells
    cellArea = Cells(c).Area;
    Ft(c,:) = (sumPixels(c, :) /cellArea);
end


%2.Calculate the time-dependent baseline F0(t )
F0 = zeros(numCells, time);
traceSmooth = zeros(numCells, time);
%smoothing
for i = 1:numCells
    traceSmooth(i,:) = smooth(Ft(i,:), 'loess');
end

for t = 1:time
    t1 = t-tau1;
    if t1 < 1
        t1 = 1;
    end
    F0(:, t) = min(traceSmooth(:,t1:t),[],2);
end

%3.Calculate the relative change of fluorescence signal R(t) from F(t) and F0(t ) =
Rt = (Ft - F0) ./F0;

%4.exponentially weighted moving average- the first tau2 is always NaN so
%correction is made below
for c = 1:numCells
    df_F0(c,:) = tsmovavg(Rt(c,:),'e',tau2, 2);
end
df_F0_mat = repmat(df_F0(:,tau2), [1, tau2-1]);

df_F0(:,1:tau2-1) = df_F0_mat;

%% new trace - neuropil

df_F0_neuropil = df_F0 - traceNeuropilSmooth;

% %2.Calculate the time-dependent baseline F0(t )
% F0 = zeros(numCells, time);
% traceSmooth = zeros(numCells, time);
% %smoothing
% for i = 1:numCells
%     traceSmooth(i,:) = smooth(Ftneuropil(i,:), 'loess');
% end
% 
% for t = 1:time
%     t1 = t-tau1;
%     if t1 < 1
%         t1 = 1;
%     end
%     F0(:, t) = min(traceSmooth(:,t1:t),[],2);
% end
% 
% %3.Calculate the relative change of fluorescence signal R(t) from F(t) and F0(t ) =
% Rt = (Ftneuropil - F0) ./F0;
% 
% %4.exponentially weighted moving average- the first tau2 is always NaN so
% %correction is made below
% df_F0_neuropil = zeros(numCells, time);
% for c = 1:numCells
%     df_F0_neuropil(c,:) = tsmovavg(Rt(c,:),'e',tau2, 2);
% end
% df_F0_mat = repmat(df_F0_neuropil(:,tau2), [1, tau2-1]);
% 
% df_F0_neuropil(:,1:tau2-1) = df_F0_mat;

%% Cross Correlation method with several spike examples
load('spikeShapes2.mat');
localMax = 10;
scale_SD = 1;
locs3=cell(numCells,1);
traceChosen = df_F0;
for i = 1:numCells
    %normxcorr2
    meanNormalSpike = meanNormalSpike - mean(meanNormalSpike);
    meanFastSpike = meanFastSpike - mean(meanFastSpike);
    meanFlatSpike = meanFlatSpike - mean(meanFlatSpike);
    twoSpikes = twoSpikes - mean(twoSpikes);
    [tracexCorrTemp{1}] = normxcorr2(meanNormalSpike, traceChosen(i,:));
    [tracexCorrTemp{2}] = normxcorr2(meanFastSpike, traceChosen(i,:));
    [tracexCorrTemp{3}] = normxcorr2(meanFlatSpike, traceChosen(i,:));
    [tracexCorrTemp{4}] = normxcorr2(twoSpikes, traceChosen(i,:));
    for t = 1:length(tracexCorrTemp)
        tracexCorr = tracexCorrTemp{t};
        thresh = mean(tracexCorr) + std(tracexCorr)*scale_SD;
        traceThresh = tracexCorr.*(tracexCorr > thresh);
        [pks,locsxCorr] = findpeaks(traceThresh);
        lagDiff = length(meanFastSpike)-1;
        %[~,I] = max(abs(locsxCorr));
        %lagDiff = lags(I);
        locsTemp = locsxCorr - lagDiff;
        locsTemp(locsTemp > time) = time;
        locsTemp(locsTemp < 1) = 1;
        
        %check threshold of original trace as well
        if ~isempty(locsTemp)
            %shift to local max
            [locs] = findLocalMax(locsTemp, traceChosen, i, time, localMax);
            %check that all peaks are over 0
            threshOrigTrace = mean(traceChosen(i,:)) + std(traceChosen(i,:))*scale_SD;
            locs(traceChosen(i, locs) < threshOrigTrace) = [];
            zeroThresh = 0;
            locs(traceChosen(i, locs) < zeroThresh) = [];
            [locs] = removeCloseSpikes(locs);
            locs3{i,1} = sort([locs3{i}, locs], 'ascend');
        end
    end
end

for l = 1:length(locs3)
    [locs3{l,1}] = removeCloseSpikes(locs3{l});
    numSpikes(l) = length(locs3{l});
end

%% Plot
figure;
subplot(4,1,1);
plot(df_fixedF0(1, :));
title('Trace for min10 Df/F0')

subplot(4,1,2);
plot(df_fixedF0Back(1, :));
title('Trace for min10 Df/F0 neuropil subtracted')

subplot(4,1,3);
plot(df_F0(1, :));
hold on
plot(locs3{1},df_F0(1,locs3{1}), 'r*');
title('New trace');
hold off

subplot(4,1,4);
plot(df_F0_neuropil(1, :));
title('New trace -neuropil');

%% save
handles.df_fixedF0 = traceChosen;
%handles.df_fixedF0WOBack=df_fixedF0WOBack;
handles.locs = locs3;
handles.numSpikes = numSpikes;

end
