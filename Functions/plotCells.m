function [] = plotCells(handles)
%Plots new cells and points


Cells = handles.Cells;
cellLocationsCells = [];
n = handles.n;
image = handles.meanImage;
activeCells = handles.activeCells;
ITimeseriesSTD = handles.ITimeseriesSTD;

binaryImage = zeros(handles.sizeImage);
binaryImageActiveCells = zeros(size(ITimeseriesSTD,1),size(ITimeseriesSTD,1));
cellLocationsCells=zeros(length(Cells),2);
for c = 1:length(Cells)
    cellLocationsCells(c,:) = Cells(c).Centroid;
    index = Cells(c).PixelIdxList;
    binaryImage(index)=c;
end
cellLocationsActiveCells=zeros(length(activeCells),2);

for c = 1:length(activeCells)
    cellLocationsActiveCells(c,:) = activeCells(c).Centroid;
end


% Plot interest points
if ~isempty(cellLocationsCells)
    axes(handles.axes1);
    imagesc(image); colormap(gray); hold on;
    title('Mean Image with segmented cells');
    axis off;
    hold on;
    
    [B,L] = bwboundaries(binaryImage,'noholes');
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1)
    end
    for i = 1:size(cellLocationsCells,1)
        hnd1=text(cellLocationsCells(i,1),cellLocationsCells(i,2),num2str(i));
        set(hnd1,'FontSize',12, 'Color', 'w')

    end
    hold off;

    axes(handles.axes2);
    imagesc(ITimeseriesSTD); colormap(gray);
    title('STD Image with active cells');
    axis off;
    hold on;
    
    plot(cellLocationsActiveCells(:,1),cellLocationsActiveCells(:,2), 'w*');
    
    for i = 1:size(cellLocationsActiveCells,1)
        hnd1=text(cellLocationsActiveCells(i,1),cellLocationsActiveCells(i,2),num2str(i));
        set(hnd1,'FontSize',12, 'Color', 'w')
    end

    hnd1=text(cellLocationsActiveCells(n,1),cellLocationsActiveCells(n,2),num2str(n));
    set(hnd1,'FontSize',12, 'Color', 'r')

    hold off;
else
    axes(handles.axes1);
    imagesc(image); colormap(gray); hold on;
    title('Mean Image with segmented cells');
    axis off;
    hold off;

    axes(handles.axes2);
    imagesc(ITimeseriesSTD); colormap(gray);
    title('STD Image with active cells');
    axis off;
    hold on;

    display('nothing to plot');

end
end