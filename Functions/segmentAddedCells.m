function [Cells] = segmentAddedCells(I, handles, cellLocations)
%This function recieves mean image detects points, and segments them
Plot = 0;
cellSize = 30;
neuropilSize = 50;
con = handles.con;
sizeImage = size(I,1);
numCellLocs = size(cellLocations,1);
eccThresh = 0.95;

sizeITs = 256;

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
    BW2 = patch > thresh;
    
    %remove small regions
    %BW2 = bwareaopen(BW,100,con);
    
    %dilate
    %BW2 = bwmorph(BW,'close',inf);  
    
    %remove all but 1 cells
    [BWLabelled, num] = bwlabel(BW2, con);
    if num > 1
        area = [];
        tempPatch = zeros([size(patch),num]);
        for j = 1:num
           tempPatch(:,:,j) = BWLabelled == j;
           area(j) = bwarea(tempPatch(:,:,j));
        end
    
    [~,index1] = max(area);
    BW2 = tempPatch(:,:,index1)> 0;
    end
    
    binaryImage(y1(n):y2(n),x1(n):x2(n)) = BW2;
    
    if Plot==1
        figure, imagesc(patch); colormap(gray);
        title('Orignal patch');
        movegui('west');
        axis off;
        figure, imagesc(BW2); colormap(gray);
        title('Segmented Image');
        movegui('east');
        axis off;        
    end
end

%closing
%binaryImage = bwmorph(binaryImage,'close',inf); 

%remove small regions
%binaryImage = bwareaopen(binaryImage,100,con);

%remove high elonaged structures using Eccentricity
Cells=regionprops(binaryImage, 'Area', 'Centroid', 'Eccentricity', 'PixelIdxList', 'BoundingBox');


%Get neuropil mask per cell
for c = 1:length(Cells)
    tempImageCell = zeros(sizeImage);
    tempImageBox= zeros(sizeImage);
    tempImageCell(Cells(c).PixelIdxList) = 1;
    Centroid = Cells(c).Centroid;
    x1 = round(Centroid(1) - round(neuropilSize/2));
    x2 = round(Centroid(1) + round(neuropilSize/2));
    y1 = round(Centroid(2) - round(neuropilSize/2));
    y2 = round(Centroid(2) + round(neuropilSize/2));

    %Ensure indeces are within the image margins
    if x1 <= 0
        x1 = 1;
    end
    if y1 <= 0
        y1 = 1;
    end
    if x2 > sizeImage
        x2 = round(sizeImage);
    end
    if y2 > sizeImage
        y2 = round(sizeImage);
    end
    
    tempImageBox(y1:y2, x1:x2) = 1;
    neuropilImage = (tempImageBox - tempImageCell)>0;
    tempNeuropilPixels=regionprops(neuropilImage, 'PixelIdxList');
    Cells(c).neuropilPixelIdxList = tempNeuropilPixels(1).PixelIdxList;
end
end