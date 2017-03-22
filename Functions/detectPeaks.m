function [handles] = detectPeaks(trace, traceNeuropil, handles)

Plot = 0;
validMask = [1 1 1 0 -1 -1 -1];
numCells = size(trace,1);
time = size(trace,2);
thresh = mean(trace(:))+std(trace(:));


%% df/d calculation
sort_trace_by_F = sort(trace, 2, 'ascend');
for x = 1:numCells
min_10(x, 1) = mean(sort_trace_by_F(x, 1:10));
end

%% (F-F0) / F0 
df_fixedF0WOBack = zeros(numCells,time);
for c=1:numCells
    for t=1:(time-1)
        df_fixedF0WOBack(c,t) = (trace(c,t+1)- min_10(c))/(min_10(c));
    end
    df_fixedF0WOBack(c,t) = (trace(c,t+1)- min_10(c))/(min_10(c));
end
df_fixedF0WOBack(:,time) = df_fixedF0WOBack(:,time-1);

%% df/d calculation using the neuropil (i.e. surrounding area) to normalize

%with smoothing
for i = 1:numCells
    traceNeuropilSmooth(i,:) = smooth(traceNeuropil(i,:), 'loess');
end
min10mat = repmat(min_10, [1, time]);
df_fixedF0 = ((trace - traceNeuropilSmooth)-min10mat)./min10mat;
%Plot neuropil traces
% figure;
% subplot(3,1,1);
% plot(1:1:handles.time, trace(17,:));
% title('Plot of raw trace');
% subplot(3,1,2);
%plot(1:1:handles.time, df_fixedF0(12,:));
%title('Plot of DF/F0 trace- with background substracted');
% subplot(3,1,3);
% plot(1:1:handles.time, df_fixedF0WOBack(17,:));
% title('Plot of DF/F0 trace- with background substracted');
% pause(1);

%% Cross Correlation method with several spike examples
load('spikeShapes2.mat');
localMax = 10;
scale_SD = 2;
locs3=cell(numCells,1);
for i = 1:numCells
    %normxcorr2
    meanNormalSpike = meanNormalSpike - mean(meanNormalSpike);
    meanFastSpike = meanFastSpike - mean(meanFastSpike);
    meanFlatSpike = meanFlatSpike - mean(meanFlatSpike);
    twoSpikes = twoSpikes - mean(twoSpikes);
    [tracexCorrTemp{1}] = normxcorr2(meanNormalSpike, df_fixedF0WOBack(i,:));
    [tracexCorrTemp{2}] = normxcorr2(meanFastSpike, df_fixedF0WOBack(i,:));
    [tracexCorrTemp{3}] = normxcorr2(meanFlatSpike, df_fixedF0WOBack(i,:));
    [tracexCorrTemp{4}] = normxcorr2(twoSpikes, df_fixedF0WOBack(i,:));
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
            [locs] = findLocalMax(locsTemp, df_fixedF0WOBack, i, time, localMax);
            %check that all peaks are over 0
            threshOrigTrace = mean(df_fixedF0WOBack(i,:)) + std(df_fixedF0WOBack(i,:))*scale_SD;
            locs(df_fixedF0WOBack(i, locs) < threshOrigTrace) = [];
            zeroThresh = 0;
            locs(df_fixedF0WOBack(i, locs) < zeroThresh) = [];
            [locs] = removeCloseSpikes(locs);
            locs3{i,1} = sort([locs3{i}, locs], 'ascend');
        end
    end
end

for l = 1:length(locs3)
    [locs3{l,1}] = removeCloseSpikes(locs3{l});
    numSpikes(l) = length(locs3{l});
end
%% save
handles.df_fixedF0 = df_fixedF0WOBack;
%handles.df_fixedF0WOBack=df_fixedF0WOBack;
handles.locs = locs3;
handles.numSpikes = numSpikes;

end