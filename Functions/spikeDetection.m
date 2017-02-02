

trace = Data.trace;
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

for i = 1:numCells
    locs1(i,:) = artifactThresh(df_fixedF0(1,:),validMask,thresh);
    numSpikes = length(locs1);
    figure;
    plot(1:1:800,df_fixedF0(i,:));
    title(['Number of peaks detected: ', num2str(numSpikes)]);
    hold on
    plot(locs1(i,:),df_fixedF0(i,locs1(i,:)), 'r*');
    pause(1);
end

%% Threshold method
for i = 1:numCells
    thresh = mean(df_fixedF0(i,:))+std(df_fixedF0(i,:));
    locs2{i} = find(df_fixedF0(i,:) > thresh);
    numSpikes = length(locs2);
    figure;
    plot(1:1:800,df_fixedF0(i,:));
    title(['Number of peaks detected: ', num2str(numSpikes)]);
    hold on
    plot(locs2{i},df_fixedF0(i,locs2{i}), 'r*');
    pause(1);
end
