function [] = plotPeaks(handles)
n = handles.n;
axes(handles.axes3);
df_fixedF0 = handles.df_fixedF0;

locs = handles.locs;
    plot(1:1:handles.time,df_fixedF0(n,:));
    title(['Number of cells detected: ', num2str(length(locs))]);
if ~isempty(locs{n})
    hold on 
    plot(locs{n},df_fixedF0(n,locs{n}), 'r*');
    hold off
end


end