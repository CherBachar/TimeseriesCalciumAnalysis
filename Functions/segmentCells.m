function [Cells] = segmentCells(I, handles, cellLocations)
%This function recieves mean image detects points, and segments them
Plot = 0;
cellSize = handles.cellSize;
con = handles.con;
sizeImage = handles.sizeImage;
numCellLocs = size(cellLocations,1);
eccThresh = 0.95;
%% Extract bouton patches
%Extract the bouton patches using the bouton locations extracted
for n = 1:numCellLocs
    x1(n) = round(cellLocations(n,1) - round(cellSize/2));
    x2(n) = round(cellLocations(n,1) + round(cellSize/2));
    y1(n) = round(cellLocations(n,2) - round(cellSize/2));
    y2(n) = round(cellLocations(n,2) + round(cellSize/2));

    %Ensure indeces are within the image margins
    if x1(n) <= 0
        x1(n) = 1;
    end
    if y1(n) <= 0
        y1(n) = 1;
    end
    if x2(n) > sizeImage
        x2(n) = round(sizeImage);
    end
    if y2(n) > sizeImage
        y2(n) = round(sizeImage);
    end
end

[cellPatch] = extractBoutonPatch(cellLocations, cellSize, sizeImage, I, 0);

%create a binary image
binaryImage = zeros(sizeImage);

%% Segmentation
for n = 1:length(cellPatch)
    patch = cellPatch{n};
    maxPixelVal = max(patch(:));
    minPixelVal = min(patch(:));
    thresh2 = minPixelVal + ((maxPixelVal - minPixelVal) * 0.5);
    meanPixelVal = mean(patch(:));
    stdPixel = std(patch(:));
    thresh = meanPixelVal + stdPixel;
    BW = patch > thresh;
    
    %remove small regions
    BW2 = bwareaopen(BW,100,con);
    
    %dilate
    BW2 = bwmorph(BW2,'close',inf);  
    
    %remove all but 2 cells
    [BWLabelled, num] = bwlabel(BW2, con);
    if num > 1
        area = [];
        tempPatch = zeros([size(patch),num]);
        for j = 1:num
           tempPatch(:,:,j) = BWLabelled == j;
           area(j) = bwarea(tempPatch(:,:,j));
        end
    
    [~,index1] = max(area);
    area(index1)=0;
    [~,index2] = max(area);
    BW2 = (tempPatch(:,:,index1) + tempPatch(:,:,index2))> 0;
    end
    
    binaryImage(y1(n):y2(n),x1(n):x2(n)) = BW2;
    
    
    if Plot == 1
        figure, imagesc(patch); colormap(gray);
        title('Orignal patch');
        movegui('west');
        axis off;
        figure, imagesc(BW2); colormap(gray);
        title('Segmented Image');
        movegui('east');
        axis off;
        pause(1);
    end
end

%closing
%binaryImage = bwmorph(binaryImage,'close',inf); 

%remove small regions
binaryImage = bwareaopen(binaryImage,100,con);

%remove high elonaged structures using Eccentricity
Cells=regionprops(binaryImage, 'Centroid', 'Eccentricity', 'PixelIdxList');

for r = 1:length(Cells)
    remEccentricity(r) = Cells(r).Eccentricity > eccThresh;
end

Cells(remEccentricity) = [];

end