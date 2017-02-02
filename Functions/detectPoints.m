function [handles] = detectPoints(I, handles)
Plot = 0;
detector = handles.detector;
shiftCentroid = handles.shiftCentroid;
distThresh = handles.distThresh;
cellSize = handles.cellSize;
con = handles.con;
sizeImage = handles.sizeImage;

[LoGImage, padSize] = LoGConv(I, Plot); %convolve image with LoG mask
LoGNorm = LoGImage / max((LoGImage(:)));

if strcmp(detector, 'SURF')
    %SURF keypoints
    SURFPoints = detectSURFFeatures(LoGNorm);
    SURFLocations = SURFPoints.Location;
    SURFLocations = SURFLocations - padSize(1);
    detectorLocations.SURF = SURFLocations;
    interestPoints = SURFLocations;
    
elseif strcmp(detector, 'harris')
    %Harris keypoints
    harrisPoints = detectHarrisFeatures(LoGNorm);
    harrisLocations = harrisPoints.Location;
    harrisLocations = harrisLocations - padSize(1);
    detectorLocations.harris = harrisLocations;
    interestPoints = harrisLocations;
    
elseif strcmp(detector, 'SIFT')
    %SIFT keypoints
    [SIFTPoints, SIFTFeatures] = vl_sift(single(LoGNorm));
    [~, indexFeatures] = sort(sum(SIFTFeatures,1), 'ascend');
    SIFTLocations = SIFTPoints(1:2,indexFeatures(1:100))';
    SIFTLocations = SIFTLocations - padSize(1);
    detectorLocations.SIFT = SIFTLocations;
    interestPoints = SIFTLocations;
end    

%check all points are positive
[row,~] = find(interestPoints < 1);
interestPoints(row,:) = [];
[row,~] = find(interestPoints > sizeImage);
interestPoints(row,:) = [];

cellLocations = interestPoints;
%% Plot interest points
if Plot==1
    figure; imagesc(I); colormap(gray); hold on;
    plot(interestPoints(:,1),interestPoints(:,2),'r+')
    title('Cell interest points');
    axis off;
end

%% Remove locations with low intensity
% 
% meanIntensity = mean(I(:));
% STDIntensity = std(double(I(:)));
% 
% for i = 1:length(interestPoints)
%     intensities(i) = I(round(interestPoints(i,1)),round(interestPoints(i,2)));
% end
% 
% [~,col] = find(intensities < (meanIntensity - STDIntensity));
% interestPoints(col,:) = [];
% 
% % Plot interest points
% if Plot==1
%     figure; imagesc(I); colormap(gray); hold on;
%     plot(interestPoints(:,1),interestPoints(:,2),'b+')
%     title('Cell interest points');
%     axis off;
% end
%% Shift centroids to local maximum

%shift interest points to their local max
%[cellLocations] = shiftCentroidsToLocalMax(cellLocations, I, shiftCentroid);

%% Discard close boutons
%remove boutons closer than distThresh. Bouton with the lower pixel
%intensity is removed

[cellLocations] = removeAllCloseBoutons(cellLocations, I, distThresh);
% Plot interest points
axes(handles.axes1);
imagesc(I); colormap(gray); hold on;
plot(cellLocations(:,1),cellLocations(:,2),'g+')
title('Cell interest points');
axis off;
hold off;
handles.cellLocations = cellLocations;
end