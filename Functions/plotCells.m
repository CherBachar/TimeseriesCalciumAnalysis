function [] = plotCells(handles)
%Plots new cells and points

Cells = handles.Cells;
binaryImage = zeros(handles.sizeImage);
cellLocations = [];
n = handles.n;

for c = 1:length(Cells)
    cellLocations(c,:) = Cells(c).Centroid;
    index = Cells(c).PixelIdxList;
    binaryImage(index)=c;
end

% Plot interest points
axes(handles.axes1);
imagesc(handles.meanImage); colormap(gray); hold on;
plot(cellLocations(:,1),cellLocations(:,2),'g+')
plot(cellLocations(n,1),cellLocations(n,2),'r+')
title('Mean Image with the cell centroids');
axis off;
hold off;

axes(handles.axes2);
imagesc(handles.meanImage); colormap(gray);
title('Mean Image with segmented cells');
axis off;
hold on;


[B,L] = bwboundaries(binaryImage,'noholes');
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
end
for i = 1:size(cellLocations,1)
    hnd1=text(cellLocations(i,1),cellLocations(i,2),num2str(i));
    set(hnd1,'FontSize',12, 'Color', 'w')

end

hnd1=text(cellLocations(n,1),cellLocations(n,2),num2str(n));
set(hnd1,'FontSize',12, 'Color', 'r')

hold off;
end