function [df_F0] = df_F0Trace(trace, handles)

%% initialize
numCells = size(trace,1);
time = size(trace,2);
thresh = mean(trace(:))+std(trace(:));
sizeImage = handles.sizeImage;
Cells = handles.activeCells;
sumPixels = handles.sumPixels;
tau1 = 10;
tau2 = 10;

%% df/f calculation from paper: In vivo two-photon imaging of sensory-evoked dendritic calcium signals in cortical neurons
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

%% save
handles.df_F0 = df_F0;

end