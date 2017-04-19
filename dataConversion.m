%% Conversion data script
%script for converting data from 'Legacy' version to the 'master' version
%of this software- saves new version in current folder with the filename

%open files for conversion- load all files for conversion
[fileName,PathName] = uigetfile('*.mat', 'MultiSelect','on');

if ~iscell(fileName)
    FileName{1} = fileName;
else 
    FileName = fileName;
end
for i = 1:length(FileName)
    load([PathName, FileName{i}]);
    
    %Get Area and put in Cells
    Cells = Data.Cells;
    sizeImage = size(Data.meanImage,1);
    binaryImage = zeros(sizeImage);
    %Cells and neuropil
    for l = 1:length(Cells)
        %Get cell mask
        tempMaskCell = zeros(sizeImage);
        index = Cells(l).PixelIdxList;
        binaryImage(index)=1;
        tempMaskCell(index)=1;
        AreaTemp = regionprops(tempMaskCell, 'Area');
        Cells(l).Area = AreaTemp.Area;
    end
    
    %Reanalyse 
    ITimeseries =Data.ITimeseries;
    numLayers = size(ITimeseries,3); 
    sizeITs = size(ITimeseries,1);

    binaryImage = zeros(sizeImage);
    binaryImageNeuropil = zeros(sizeImage);

    %Cells and neuropil
    for l = 1:length(Cells)
        %Get cell mask
        tempMaskCell = zeros(sizeImage);
        index = Cells(l).PixelIdxList;
        binaryImage(index)=1;
        tempMaskCell(index)=1;
        maskResized(:,:,l) = imresize(tempMaskCell, [sizeITs sizeITs])>0;
        %Get neuropil mask
        tempMaskNeuropil = zeros(sizeImage);
        index = Cells(l).neuropilPixelIdxList;
        binaryImageNeuropil(index)=1;
        tempMaskNeuropil(index)=1;
        tempMaskNeuropil = (tempMaskNeuropil - tempMaskCell) >0;
        maskNeuropilResized(:,:,l) = imresize(tempMaskNeuropil, [sizeITs sizeITs])>0;
    end

    binaryImageNeuropil = (binaryImageNeuropil - binaryImage) >0;
    binaryImageResized = imresize(binaryImage, [sizeITs sizeITs]);
    binaryImageNeuropilResized = imresize(binaryImageNeuropil, [sizeITs sizeITs]);

    [B,L] = bwboundaries(binaryImageResized,'noholes');
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
    end
    hold off;


    disp('in progress');
    for c = 1:length(Cells)
        maskTempCell = logical(maskResized(:,:,c)); 
        maskTempNeuropil = logical(maskNeuropilResized(:,:,c)); 

        for time = 1:numLayers
            tempIT1 = ITimeseries(:,:,time);
            trace(c,time) = mean(mean(tempIT1(maskTempCell)));
            traceNeuropil(c,time) = mean(mean(tempIT1(maskTempNeuropil)));
        end    
        disp(['Cell: ', num2str(c)]);

    end

    handles.trace = trace;
    handles.time = numLayers;
    handles.ITimeseries=ITimeseries;

    % Detect peaks on traces
    [handles] = detectPeaks(trace, traceNeuropil, handles);
    df_fixedF0 = handles.df_fixedF0;

    %Find synchrony 
    [handles] = findCaSynchrony(df_fixedF0, handles);
    disp('Finished analysis');

    %save
    Data.Cells = Cells;
    Data.ITimeseries = handles.ITimeseries;
    Data.trace = handles.trace;
    Data.df_fixedF0 = handles.df_fixedF0;
    Data.locs = handles.locs;
    Data.numSpikes = handles.numSpikes;
    Data.synchrony = handles.synch; 
    
    save([FileName{i}], 'Data');
    
    clearvars -except FileName PathName
    close all
    
end