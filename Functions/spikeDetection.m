

trace = Data.trace;
validMask = [1 1 1 0 -1 -1 -1];
numCells = size(trace,1);
time = size(trace,2);
thresh = mean(trace(:))+std(trace(:));


%% df/d calculation
df_fixedF0 = zeros(time,numCells);
min_10 = zeros(1,numCells);
trace2 = trace';
for r=1:time
for c=(1:numCells-1)
df_fixedF0(r,c) = (trace2(r,c+1)- min_10(1,c))/(min_10(1,c));
end
end

%% using https://github.com/LeventhalLab/EphysToolbox/tree/master/SpikeySpike

% for i = 1:numCells
%     locs(i,:) = artifactThresh(trace(1,:),validMask,thresh);
%     numSpikes = length(locs);
%     figure;
%     plot(1:1:800,trace(i,:));
%     title(['Number of cells detected: ', num2str(numSpikes)]);
%     hold on
%     plot(locs(i,:),trace(i,locs(i,:)), 'r*');
%     pause(1);
% end

%% Threshold method
trace3 = df_fixedF0';
for i = 1:numCells
    %thresh = mean(trace3(i,:))+std(trace3(i,:));
    %locs{i} = find(trace3(i,:) > thresh);
    %numSpikes = length(locs);
    figure;
    plot(1:1:800,trace3(i,:));
    %title(['Number of cells detected: ', num2str(numSpikes)]);
    hold on
    %plot(locs{i},trace3(i,locs{i}), 'r*');
    pause(1);
end
