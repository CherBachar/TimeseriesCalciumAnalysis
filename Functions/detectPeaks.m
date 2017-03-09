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
% plot(1:1:handles.time, df_fixedF0(17,:));
% title('Plot of DF/F0 trace- with background substracted');
% subplot(3,1,3);
% plot(1:1:handles.time, df_fixedF0WOBack(17,:));
% title('Plot of DF/F0 trace- with background substracted');
% 
% pause(1);

%% DF/F0 using kalman filter

%% using https://github.com/LeventhalLab/EphysToolbox/tree/master/SpikeySpike
%http://gaidi.ca/weblog/extracting-spikes-from-neural-electrophysiology-in-matlab
% for i = 1:numCells
%     locs1(i,:) = artifactThresh(df_fixedF0(1,:),validMask,thresh);
%     numSpikes = length(locs1);
%     if Plot == 1
%         figure;
%         plot(1:1:time,df_fixedF0(i,:));
%         title(['Number of peaks detected: ', num2str(numSpikes)]);
%         hold on
%         plot(locs1(i,:),df_fixedF0(i,locs1(i,:)), 'r*');
%         pause(1);
%     end
% end

%% Threshold method
% for i = 1:numCells
%     thresh = mean(df_fixedF0(i,:))+std(df_fixedF0(i,:));
%     locs2{i} = find(df_fixedF0(i,:) > thresh);
%     numSpikes(i) = length(locs2);
%     if Plot == 1
%         figure;
%         plot(1:1:time,df_fixedF0(i,:));
%         title(['Number of peaks detected: ', num2str(numSpikes(i))]);
%         hold on
%         plot(locs2{i},df_fixedF0(i,locs2{i}), 'r*');
%         pause(1);
%     end
% end

%% Cross-correlation method
% load('CrossCorr.mat');
% localMax = time / 80;
% for i = 1:numCells
%     [tracexCorr,lags] = xcorr(df_fixedF0(i,:),CrossCorrMask,'none');
%     thresh = mean(tracexCorr) + std(tracexCorr);
%     traceThresh = tracexCorr.*(tracexCorr > thresh);
%     [pks,locsxCorr] = findpeaks(traceThresh);
%     
%     [~,I] = max(abs(locsxCorr));
%     lagDiff = lags(I);
%     %s1al = df_fixedF0(-lagDiff+1:end);
%     locsTemp1 = locsxCorr + lagDiff;
%     locsTemp1(locsTemp1 > time) = time;
%     locsTemp = [];
%     for l = 1:length(locsTemp1)
%         x1 = locsTemp1(l)-localMax;
%         x2 = locsTemp1(l)+localMax;
%         if (locsTemp1(l)+localMax) > time
%             x2 = time;
%         end
%         if (locsTemp1(l)-localMax) < 1
%             x1 = 1;
%         end
%         [~,ind] = max(df_fixedF0(i,x1:x2));
%         locsTemp(l) = locsTemp1(l) - localMax - 1 + ind;
%         
%         %make sure it's between 1-time
%         if locsTemp(l) < 1
%             locsTemp(l) = 1;
%         end
%         if locsTemp(l) > time
%             locsTemp(l) = time;
%         end        
%     end
%     locs3{i} = locsTemp;
%     numSpikes(i) = length(locsTemp);
% end
%% Cross Correlation method with several spike examples
load('spikeShapes.mat');
localMax = 10;
scale_SD = 2;
locs3=cell(numCells,1);
for i = 1:numCells
    [tracexCorrTemp{1},lags] = xcorr(df_fixedF0(i,:),meanNormalSpike,'none');
    [tracexCorrTemp{2},lags] = xcorr(df_fixedF0(i,:),meanFastSpike,'none');
    [tracexCorrTemp{3},lags] = xcorr(df_fixedF0(i,:),meanFlatSpike,'none');
    for t = 1:length(tracexCorrTemp)
        tracexCorr = tracexCorrTemp{t};
        thresh = mean(tracexCorr) + std(tracexCorr)*scale_SD;
        traceThresh = tracexCorr.*(tracexCorr > thresh);
        [pks,locsxCorr] = findpeaks(traceThresh);

        [~,I] = max(abs(locsxCorr));
        lagDiff = lags(I);
        locsTemp = locsxCorr + lagDiff;
        locsTemp(locsTemp > time) = time;

        [locs] = findLocalMax(locsTemp, df_fixedF0, i, time, localMax);
        
        locs3{i,1} = sort([locs3{i}, locs], 'ascend');
    end
end

for l = 1:length(locs3)
    [locs3{l,1}] = removeCloseSpikes(locs3{l});
    numSpikes(l) = length(locs3{l});
    
end
%% save
handles.df_fixedF0 = df_fixedF0;
handles.df_fixedF0WOBack=df_fixedF0WOBack;
handles.locs = locs3;
handles.numSpikes = numSpikes;

end