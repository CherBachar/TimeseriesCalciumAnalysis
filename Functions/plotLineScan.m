function [] = plotLineScan(handles, activeCellsMask)
%Plot time-series line scan of neuropil and cells

ITimeseries = handles.ITimeseries;
activeCells = handles.activeCells;
neuropilMask = handles.neuropilMask;
numCells = length(activeCells);
time = size(ITimeseries,3);
cellScan = {};
for c = 1:numCells
    tempCellMask = activeCellsMask(:,:,c);
    cellScanTemp=[];
    for t = 1:time
        temp = ITimeseries(:,:,t);
        cellScanTemp(:,t)=temp(logical(tempCellMask));
        neuropilScan(:,t) = temp(neuropilMask);
        cellScan{c,1} = cellScanTemp;
    end
    

end

%sample
% n_neuropil = size(neuropilScan,1);
% neuropil_index = round(linspace(1,n_neuropil,1000));
% neuropilScanSample = neuropilScan(neuropil_index,:);

figure;
subplot(numCells+1,1,1);
imagesc(mat2gray(neuropilScan));
axis off;
title('Line scan plot of the neuropil','Fontsize', 7);

for c = 1:numCells
    subplot(numCells+1,1,1+c);
    plotCell = cellScan{c};
    imagesc(mat2gray(plotCell));
    title(['Cell ', num2str(c)], 'Fontsize', 7);
    axis off;
end

hp4 = get(subplot(numCells+1,1,1+c),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.05  hp4(2)+hp4(3)-0.05])
    
% matrixNeuropilCells = neuropilScan;
% for c = 1:numCells
%     matrixNeuropilCells = [matrixNeuropilCells; cellScan{c}]; 
% end
% figure;imagesc(matrixNeuropilCells);colorbar;
end