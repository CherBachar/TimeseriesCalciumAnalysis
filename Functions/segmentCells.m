function [handles] = segmentCells(I, handles)
%This function recieves mean image detects points, and segments them
Plot = 0;
%detector = handles.detector;
%shiftCentroid = handles.shiftCentroid;
%distThresh = handles.distThresh;
cellSize = handles.cellSize;
con = handles.con;
sizeImage = handles.sizeImage;
cellLocations = handles.cellLocations;

% [LoGImage, padSize] = LoGConv(I, Plot); %convolve image with LoG mask
% LoGNorm = LoGImage / max((LoGImage(:)));
% 
% if strcmp(detector, 'SURF')
%     %SURF keypoints
%     SURFPoints = detectSURFFeatures(LoGNorm);
%     SURFLocations = SURFPoints.Location;
%     SURFLocations = SURFLocations - padSize(1);
%     detectorLocations.SURF = SURFLocations;
%     interestPoints = SURFLocations;
%     
% elseif strcmp(detector, 'harris')
%     %Harris keypoints
%     harrisPoints = detectHarrisFeatures(LoGNorm);
%     harrisLocations = harrisPoints.Location;
%     harrisLocations = harrisLocations - padSize(1);
%     detectorLocations.harris = harrisLocations;
%     interestPoints = harrisLocations;
%     
% elseif strcmp(detector, 'SIFT')
%     %SIFT keypoints
%     [SIFTPoints, SIFTFeatures] = vl_sift(single(LoGNorm));
%     [~, indexFeatures] = sort(sum(SIFTFeatures,1), 'ascend');
%     SIFTLocations = SIFTPoints(1:2,indexFeatures(1:100))';
%     SIFTLocations = SIFTLocations - padSize(1);
%     detectorLocations.SIFT = SIFTLocations;
%     interestPoints = SIFTLocations;
% end    
% 
% %check all points are positive
% [row,~] = find(interestPoints < 1);
% interestPoints(row,:) = [];
% [row,~] = find(interestPoints > sizeImage);
% interestPoints(row,:) = [];
% 
% cellLocations = interestPoints;
% %% Plot interest points
% if Plot==1
%     figure; imagesc(I); colormap(gray); hold on;
%     plot(interestPoints(:,1),interestPoints(:,2),'r+')
%     title('Cell interest points');
%     axis off;
% end
% 
% %% Remove locations with low intensity
% % 
% % meanIntensity = mean(I(:));
% % STDIntensity = std(double(I(:)));
% % 
% % for i = 1:length(interestPoints)
% %     intensities(i) = I(round(interestPoints(i,1)),round(interestPoints(i,2)));
% % end
% % 
% % [~,col] = find(intensities < (meanIntensity - STDIntensity));
% % interestPoints(col,:) = [];
% % 
% % % Plot interest points
% % if Plot==1
% %     figure; imagesc(I); colormap(gray); hold on;
% %     plot(interestPoints(:,1),interestPoints(:,2),'b+')
% %     title('Cell interest points');
% %     axis off;
% % end
% %% Shift centroids to local maximum
% 
% %shift interest points to their local max
% %[cellLocations] = shiftCentroidsToLocalMax(cellLocations, I, shiftCentroid);
% 
% %% Discard close boutons
% %remove boutons closer than distThresh. Bouton with the lower pixel
% %intensity is removed
% 
% [cellLocations] = removeAllCloseBoutons(cellLocations, I, distThresh);
% % Plot interest points
% axes(handles.axes1);
% imagesc(I); colormap(gray); hold on;
% plot(cellLocations(:,1),cellLocations(:,2),'g+')
% title('Cell interest points');
% axis off;
% hold off;
%% Extract bouton patches
%Extract the bouton patches using the bouton locations extracted
for n = 1:length(cellLocations)
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
    %cc = bwconncomp(BW); %not needed
    
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

binaryImage = bwmorph(binaryImage,'close',1); 
[binaryImageLabelled, numCells] = bwlabel(binaryImage, con);

%% Plot
%IResized = imresize(I, [sizeITs sizeITs]);
%binaryImageResized = imresize(binaryImage, [sizeITs sizeITs]);
%[binaryImageLabelledResized, numCells] = bwlabel(binaryImageResized, con);

%handles.binaryImageLabelledResized = binaryImageLabelledResized;
handles.numCells = numCells;
handles.binaryImage = binaryImage;

axes(handles.axes2);
imagesc(I); colormap(gray);
title('Mean Image with segmented cells');
axis off;
hold on;

[B,L] = bwboundaries(binaryImageLabelled,'noholes');
%imshow(label2rgb(L, @jet, [.5 .5 .5]))
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end
hold off;
end