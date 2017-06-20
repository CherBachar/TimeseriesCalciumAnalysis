function [handles] = activityIntegral(handles)

%% Calculate 
% thresh = 0.5;
trace = handles.df_F0;
%thresh = mean(trace(:))+2*std(trace(:));
% traceThresh = trace.*(trace > thresh);
numCells = size(trace,1);
activityIntegral = zeros(size(trace,1),1);
for c = 1:numCells
    %calculate threshold for each cell trace
    thresh = mean(trace(c,:))+2*std(trace(c,:));
    traceThresh(c,:) = trace(c,:).*(trace(c,:)>thresh);
    activityIntegral(c) = trapz(traceThresh(c,:));
end
%% calculate synchrony between traces

% corrMatrix = zeros(numCells,numCells);
% for i = 1:numCells
%     for j = 1:numCells
%         R = corrcoef(traceThresh(i,:),traceThresh(j,:));
%         corrMatrix(i,j) = R(1,2);
%     end
% end

corrMatrix = corrcoef(traceThresh');
%% save
handles.activity = activityIntegral;
handles.corrMatrix = corrMatrix;

end