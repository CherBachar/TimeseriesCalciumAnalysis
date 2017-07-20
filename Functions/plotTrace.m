function [] = plotTrace(handles)
n = handles.n;
axes(handles.axes3);
df_F0 = handles.df_F0;
if (sum(strcmp(fieldnames(handles), 'activityThreshold')) == 1)
    thresh = handles.activityThreshold;
    thresh_n = ones(size(df_F0,2),1) .* thresh(n);
    plot(1:1:handles.time, df_F0(n,:), '-b', 1:1:handles.time, thresh_n, '--r');
else
    plot(1:1:handles.time, df_F0(n,:), '-b');
end
title(['Number of active cells detected: ', num2str(length(handles.activeCells))]);

plotCells(handles);
end