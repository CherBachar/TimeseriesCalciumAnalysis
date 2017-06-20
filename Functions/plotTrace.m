function [] = plotTrace(handles)
n = handles.n;
axes(handles.axes3);
df_F0 = handles.df_F0;

plot(1:1:handles.time,df_F0(n,:));
title(['Number of active cells detected: ', num2str(length(handles.activeCells))]);

plotCells(handles);
end