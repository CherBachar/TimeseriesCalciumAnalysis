function [handles] = findActiveCells(handles)
%Get active cells
maskResized = handles.maskResized;
ITimeseriesSTD = handles.ITimeseriesSTD;

Cells = handles.Cells;
numCells=length(Cells);
imageSize=size(ITimeseriesSTD,1);
activeCellsBinaryImage = zeros(imageSize,imageSize);
cellSize = 50;
resizeFactor = handles.sizeImage/imageSize;
%% 
a=1;
for c = 1:numCells
    
    temp=regionprops(logical(maskResized(:,:,c)), 'Centroid', 'Area', 'BoundingBox');
    num = length(temp);
    if num > 1
        area = [];
        for j = 1:num
           area(j) = temp(j).Area;
        end
    
    [~,index] = max(area);
    tempCell = temp(index);
    end
    
    centroid = tempCell.Centroid;
    x1 = round(centroid(1) - round(cellSize/2));
    x2 = round(centroid(1) + round(cellSize/2));
    y1 = round(centroid(2) - round(cellSize/2));
    y2 = round(centroid(2) + round(cellSize/2));
   
    %Check if out of boundary
    if x2 > imageSize
        x2 = imageSize;
    end
    if x1 < 1
        x1 = 1;
    end
    if y2 > imageSize
        y2 = imageSize;
    end
    if y1 < 1
        y1 = 1;
    end        
    
    patch = ITimeseriesSTD(y1:y2,x1:x2);
    
    %Hough transfrom for elipses
%     bestFits = ellipseDetection(patch);

%     %plot fit
%      figure;
%      imagesc(patch);
%     %ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 
    
    %Hough transform for circles
    [centers,radii, metric] = imfindcircles(patch,[1 10]);

    %ellipse(bestFits(:,3),bestFits(:,4),bestFits(:,5)*pi/180,bestFits(:,1),bestFits(:,2),'k');
    if ~isempty(metric)
        if (metric(1) > 0.25) && (radii(1) > 3)
            x_center = centers(1,1)+ x1;
            y_center = centers(1,2)+ y1;
            radius = radii(1)+2;
            width = Cells(c).BoundingBox(3)/resizeFactor/2;
            height = Cells(c).BoundingBox(4)/resizeFactor/2;

            x1 = round(x_center - radius); x2 = round(x_center + radius);
            y1 = round(y_center - radius); y2 = round(y_center + radius);

            %Check if out of boundary
            if x2 > imageSize
                x2 = imageSize;
            end
            if x1 < 1
                x1 = 1;
            end
            if y2 > imageSize
                y2 = imageSize;
            end
            if y1 < 1
                y1 = 1;
            end        

            activeCellsBinaryImage(y1:y2,x1:x2) = 1; 
            a=a+1;
        end
    end
end

activeCells = regionprops(activeCellsBinaryImage>0, 'Area', 'Centroid', 'Eccentricity', 'PixelIdxList', 'BoundingBox');
%% save
handles.activeCells = activeCells;
handles.percentActive = length(activeCells)/numCells*100;
end

