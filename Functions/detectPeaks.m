function [handles] = detectPeaks(trace, handles)

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

df_fixedF0 = zeros(numCells,time);
for c=1:numCells
    for t=1:(time-1)
        df_fixedF0(c,t) = (trace(c,t+1)- min_10(c))/(min_10(c));
    end
    df_fixedF0(c,t) = (trace(c,t+1)- min_10(c))/(min_10(c));
end
df_fixedF0(:,time) = df_fixedF0(:,time-1);

%% using https://github.com/LeventhalLab/EphysToolbox/tree/master/SpikeySpike
%http://gaidi.ca/weblog/extracting-spikes-from-neural-electrophysiology-in-matlab
for i = 1:numCells
    locs1(i,:) = artifactThresh(df_fixedF0(1,:),validMask,thresh);
    numSpikes = length(locs1);
    if Plot == 1
        figure;
        plot(1:1:time,df_fixedF0(i,:));
        title(['Number of peaks detected: ', num2str(numSpikes)]);
        hold on
        plot(locs1(i,:),df_fixedF0(i,locs1(i,:)), 'r*');
        pause(1);
    end
end

%% Threshold method
for i = 1:numCells
    thresh = mean(df_fixedF0(i,:))+std(df_fixedF0(i,:));
    locs2{i} = find(df_fixedF0(i,:) > thresh);
    numSpikes(i) = length(locs2);
    if Plot == 1
        figure;
        plot(1:1:time,df_fixedF0(i,:));
        title(['Number of peaks detected: ', num2str(numSpikes(i))]);
        hold on
        plot(locs2{i},df_fixedF0(i,locs2{i}), 'r*');
        pause(1);
    end
end

%% Cross-correlation method
load('CrossCorr.mat');
localMax = 10;
for i = 1:numCells
    [tracexCorr,lags] = xcorr(df_fixedF0(i,:),CrossCorrMask,'none');
    thresh = mean(tracexCorr) + std(tracexCorr)/2;
    traceThresh = tracexCorr.*(tracexCorr > thresh);
    [pks,locsxCorr] = findpeaks(traceThresh);
    
    [~,I] = max(abs(locsxCorr));
    lagDiff = lags(I);
    %s1al = df_fixedF0(-lagDiff+1:end);
    locsTemp1 = locsxCorr + lagDiff;
    locsTemp1(locsTemp1 > time) = time;
    locsTemp = [];
    for l = 1:length(locsTemp1)
        x1 = locsTemp1(l)-localMax;
        x2 = locsTemp1(l)+localMax;
        if (locsTemp1(l)+localMax) > time
            x2 = time;
        end
        if (locsTemp1(l)-localMax) < 1
            x1 = 1;
        end
        [~,ind] = max(df_fixedF0(i,x1:x2));
        locsTemp(l) = locsTemp1(l) - localMax - 1 + ind;
    end
    locs3{i} = locsTemp;
    numSpikes(i) = length(locsTemp);
end
%% save
handles.df_fixedF0 = df_fixedF0;
handles.locs = locs3;
handles.numSpikes = numSpikes;

end