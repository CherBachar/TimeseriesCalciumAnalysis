    function varargout = CalciumTimeSeriesAnalysisGUI(varargin)
% CALCIUMTIMESERIESANALYSISGUI MATLAB code for CalciumTimeSeriesAnalysisGUI.fig
%      CALCIUMTIMESERIESANALYSISGUI, by itself, creates a new CALCIUMTIMESERIESANALYSISGUI or raises the existing
%      singleton*.
%
%      H = CALCIUMTIMESERIESANALYSISGUI returns the handle to a new CALCIUMTIMESERIESANALYSISGUI or the handle to
%      the existing singleton*.
%
%      CALCIUMTIMESERIESANALYSISGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALCIUMTIMESERIESANALYSISGUI.M with the given input arguments.
%
%      CALCIUMTIMESERIESANALYSISGUI('Property','Value',...) creates a new CALCIUMTIMESERIESANALYSISGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CalciumTimeSeriesAnalysisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CalciumTimeSeriesAnalysisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CalciumTimeSeriesAnalysisGUI

% Last Modified by GUIDE v2.5 16-Feb-2017 14:26:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CalciumTimeSeriesAnalysisGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CalciumTimeSeriesAnalysisGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CalciumTimeSeriesAnalysisGUI is made visible.
function CalciumTimeSeriesAnalysisGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CalciumTimeSeriesAnalysisGUI (see VARARGIN)

% Choose default command line output for CalciumTimeSeriesAnalysisGUI
handles.output = hObject;
handles.detector = 'SURF';
handles.shiftCentroid = 1;
handles.distThresh = 10;
handles.cellSize = 50;
handles.con = 4;
handles.n = 1;
handles.load = 0;

axes(handles.axes1);
axis off
axes(handles.axes2);
axis off

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CalciumTimeSeriesAnalysisGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CalciumTimeSeriesAnalysisGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('GUI started');
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%uiwait(msgbox('Load a mean image (.tif)'));
[filename, filepath] = uigetfile('*.tif*');
handles.filename = filename;
handles.meanImage_fullfilepath = [filepath,filename];
I = double(imread(handles.meanImage_fullfilepath));
axes(handles.axes1);
imagesc(I); colormap(gray); title('original image'); axis off;
axes(handles.axes2);
imagesc(I); colormap(gray); title('original image'); axis off;

handles.sizeImage = size(I,1);
handles.meanImage = I;
handles.load = 0;

set(handles.imageName, 'String', ['File name: ', filename]);

guidata(hObject,handles);


% --- Executes on button press in SegmentCells.
function SegmentCells_Callback(hObject, eventdata, handles)
% hObject    handle to SegmentCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[cellLocations] = detectPoints(handles.meanImage, handles);
[Cells] = segmentCells(handles.meanImage, handles, cellLocations);
handles.Cells = Cells;

plotCells(handles);

guidata(hObject,handles);


% --- Executes on button press in Analyse.
function Analyse_Callback(hObject, eventdata, handles)
% hObject    handle to Analyse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Cells = handles.Cells;
inactive = 0;

if handles.load == 0
    uiwait(msgbox('Load the corresponding time-series tif stack'));
    [filename, filepath] = uigetfile('*.tif*');
    handles.imageStack_fullfilepath = [filepath,filename];
    imgInfo = imfinfo(handles.imageStack_fullfilepath);
    numLayers = length(imgInfo);

    for time = 1:numLayers
        ITimeseries(:,:,time) = double(imread(handles.imageStack_fullfilepath, time));
    end
else
    ITimeseries =handles.ITimeseries;
    numLayers = size(ITimeseries,3);
end
sizeITs = size(ITimeseries,1);

binaryImage = zeros(handles.sizeImage);
binaryImageNeuropil = zeros(handles.sizeImage);

%Cells and neuropil
for l = 1:length(Cells)
    %Get cell mask
    tempMaskCell = zeros(handles.sizeImage);
    index = Cells(l).PixelIdxList;
    binaryImage(index)=1;
    tempMaskCell(index)=1;
    maskResized(:,:,l) = imresize(tempMaskCell, [sizeITs sizeITs])>0;
    AreaTemp = regionprops(tempMaskCell, 'Area');
    cellArea(l) = AreaTemp.Area;
    %Get neuropil mask
    tempMaskNeuropil = zeros(handles.sizeImage);
    index = Cells(l).neuropilPixelIdxList;
    binaryImageNeuropil(index)=1;
    tempMaskNeuropil(index)=1;
    tempMaskNeuropil = (tempMaskNeuropil - tempMaskCell) >0;
    maskNeuropilResized(:,:,l) = imresize(tempMaskNeuropil, [sizeITs sizeITs])>0;
    %Plot neuropil binary image
%     figure;imagesc(maskNeuropilResized(:,:,l));colormap(gray);
%     title('Neuropil image');
%     axis off;
end

%save area
handles.cellArea = cellArea;

binaryImageNeuropil = (binaryImageNeuropil - binaryImage) >0;
binaryImageResized = imresize(binaryImage, [sizeITs sizeITs]);
binaryImageNeuropilResized = imresize(binaryImageNeuropil, [sizeITs sizeITs]);

%Plot neuropil binary image
figure;imagesc(binaryImageResized);colormap(gray);
title('Cell image');
axis off;
%Plot cell binary image
figure;imagesc(binaryImageNeuropilResized);colormap(gray);
title('Neuropil image');
axis off;


%Plot segmented cells on timeseries image
figure;imagesc(mean(ITimeseries,3));colormap(gray);
title('Segmentation on time series mean image');
axis off;
hold on;

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
    %mask3D = repmat(mask,[1,1,numLayers]); %Less efficient
    %maskOnly = mask3D .* ITimeseries;
    
%     figure;imagesc(maskTempNeuropil); colormap(gray);
%     movegui('west');
%     figure;imagesc(maskTempCell); colormap(gray);
%     movegui('east');
    
    for time = 1:numLayers
        tempIT1 = ITimeseries(:,:,time);
        trace(c,time) = mean(mean(tempIT1(maskTempCell)));
        traceNeuropil(c,time) = mean(mean(tempIT1(maskTempNeuropil)));
    end    
    disp(['Cell: ', num2str(c)]);
    
%     %Plot neuropil traces
%     subplot(2,1,1);
%     plot(1:1:handles.time, traceNeuropil(n,:));
%     title('Plot of Neuropil trace');
%     subplot(2,1,2);
%     plot(1:1:handles.time, trace(n,:));
%     title('Plot of segmented cell');
%     pause(1);
end

%Take into account only inactive portions of the neurpil 
% = anything less the 2*STD of the cell
if inactive == 1
    stdCell = std(trace,0, 2);
    for c = 1:length(Cells)
        maskTempNeuropil = logical(maskNeuropilResized(:,:,c)); 
        for time = 1:numLayers
            tempIT1 = ITimeseries(:,:,time);
            tempNeuropil = tempIT1(maskTempNeuropil);
            tempNeuropil = tempNeuropil(tempIT1(maskTempNeuropil) <= 2*stdCell(c));
            traceNeuropilInactive(c,time) = mean(mean(tempNeuropil));

        end
        display(['Neuropil: ', num2str(c)]);
    end
end
handles.trace = trace;
handles.time = numLayers;
handles.ITimeseries=ITimeseries;

% Detect peaks on traces
[handles] = detectPeaks(trace, traceNeuropil, handles);
axes(handles.axes3);
df_fixedF0 = handles.df_fixedF0;
df_fixedF0WOBack = handles.df_fixedF0WOBack;
locs = handles.locs;
plot(1:1:handles.time,df_fixedF0(1,:));
title(['Number of cells detected: ', num2str(length(Cells))]);
hold on
plot(locs{1},df_fixedF0(1,locs{1}), 'r*');
hold off

handles.n = 1;
figure;
subplot(3,1,1);
plot(1:1:handles.time, trace(handles.n,:));
title('Plot of raw trace');
subplot(3,1,2);
plot(1:1:handles.time, df_fixedF0(handles.n,:));
title('Plot of DF/F0 trace- with background substracted');
subplot(3,1,3);
plot(1:1:handles.time, df_fixedF0WOBack(handles.n,:));
title('Plot of DF/F0 trace- without background substracted');

%Find synchrony 
[handles] = findCaSynchrony(df_fixedF0, handles);

disp('Finished analysis');

guidata(hObject,handles);

% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Data.filename = handles.filename;
Data.meanImage = handles.meanImage;
Data.Cells = handles.Cells;
Data.ITimeseries = handles.ITimeseries;
Data.trace = handles.trace;
Data.df_fixedF0 = handles.df_fixedF0;
Data.locs = handles.locs;
Data.numSpikes = handles.numSpikes;
Data.synchrony = handles.synch;
Data.cellArea = handles.cellArea;
save([handles.filename, '.mat'], 'Data');



function CellNum_Callback(hObject, eventdata, handles)
% hObject    handle to CellNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CellNum as text
%        str2double(get(hObject,'String')) returns contents of CellNum as a double
n = str2double(get(hObject,'String'));
axes(handles.axes3);
df_fixedF0 = handles.df_fixedF0;
trace = handles.trace;
df_fixedF0WOBack = handles.df_fixedF0WOBack;

locs = handles.locs;
    plot(1:1:handles.time,df_fixedF0(n,:));
    title(['Number of cells detected: ', num2str(length(locs))]);
if ~isempty(locs{n})
    hold on 
    plot(locs{n},df_fixedF0(n,locs{n}), 'r*');
    hold off
end
handles.n = n;

figure;
subplot(3,1,1);
plot(1:1:handles.time, trace(handles.n,:));
title('Plot of raw trace');
subplot(3,1,2);
plot(1:1:handles.time, df_fixedF0(handles.n,:));
title('Plot of DF/F0 trace- with background substracted');
subplot(3,1,3);
plot(1:1:handles.time, df_fixedF0WOBack(handles.n,:));
title('Plot of DF/F0 trace- without background substracted');

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function CellNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Contrast.
function Contrast_Callback(hObject, eventdata, handles)
% hObject    handle to Contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imcontrast(gca);


% --- Executes on button press in AddPeaks.
function AddPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to AddPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'locs')
    locs = handles.locs;
    df_fixedF0 = handles.df_fixedF0;
    cell = handles.n;
    currLocs = locs{cell};
    [xAdd,~] = ginput;
    
    [xAdd] = localPeakDetector(round(xAdd), handles);
    
    currLocs = [currLocs, round(xAdd)];
    locs{cell,1} = currLocs;
    handles.locs = locs;
    handles.numSpikes(cell) = length(currLocs);

    %Plot new pooints
    axes(handles.axes3);
    plot(1:1:handles.time,df_fixedF0(cell,:));
    title(['Number of cells detected: ', num2str(length(locs))]);
    hold on
    plot(locs{cell},df_fixedF0(cell,locs{cell}), 'r*');
    hold off

    guidata(hObject,handles);
end

% --- Executes on button press in RemovePeaks.
function RemovePeaks_Callback(hObject, eventdata, handles)
% hObject    handle to RemovePeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'locs')
    locs = handles.locs;
    df_fixedF0 = handles.df_fixedF0;
    cell = handles.n;
    currLocs = locs{cell};
    [yRemove,~] = ginput;
    for i = 1:length(yRemove)
        [~, tempRemove(i)] = min(abs(currLocs - yRemove(i)));
    end
    currLocs(tempRemove) = [];
    locs{cell,1} = currLocs;
    handles.locs = locs;
    handles.numSpikes(cell) = length(currLocs);
    
    %Plot new pooints
    axes(handles.axes3);
    plot(1:1:handles.time,df_fixedF0(cell,:));
    title(['Number of cells detected: ', num2str(length(locs))]);
    hold on
    plot(locs{cell},df_fixedF0(cell,locs{cell}), 'r*');
    hold off

    guidata(hObject,handles);
end


% --- Executes on button press in AddCells.
function AddCells_Callback(hObject, eventdata, handles)
% hObject    handle to AddCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Cells = handles.Cells;
%Get new cells
[xAdd,yAdd] = ginput;
locationsAdd = [xAdd,yAdd];
if ~isempty(locationsAdd)
    locationsAdd = addPoints(locationsAdd, handles.meanImage);
    
    [newCells] = segmentCells(handles.meanImage, handles, locationsAdd);
    Cells = [Cells; newCells];
    handles.Cells = Cells;
    plotCells(handles);
end
guidata(hObject,handles);


% --- Executes on button press in RemoveCells.
function RemoveCells_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Read cell locations
Cells = handles.Cells;
for c = 1:length(Cells)
    cellLocations(c,:) = Cells(c).Centroid;
end

%Get cells to remove
[xRemove,yRemove] = ginput;
if ~isempty(xRemove)
    [~, removedCells] = removeCells(cellLocations, xRemove, yRemove);
    
    Cells(removedCells) = [];
    handles.Cells = Cells;

    plotCells(handles);
end

guidata(hObject,handles);


% --- Executes on button press in DrawBox.
function DrawBox_Callback(hObject, eventdata, handles)
% hObject    handle to DrawBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image = handles.meanImage;
sizeImage = handles.sizeImage;
box = getrect(); %box for neuropil
x1 = round(box(1));
y1 = round(box(2));
x2 = round(box(3) + x1);
y2 = round(box(4) + y1);
if x1 < 1
    x1 = 1;
end
if y1 < 1 
    y1 = 1;
end
if x2 > sizeImage
    x2 = round(sizeImage);
end
if y2 > sizeImage
    y2 = round(sizeImage);
end

neuropil = mean(mean(image(y1:y2, x1:x2)));
disp('done');
    


% --- Executes on button press in loadData.
function loadData_Callback(hObject, eventdata, handles)
% hObject    handle to loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, filepath] = uigetfile('*.mat*');
load([filepath, filename]);
set(handles.imageName, 'String', ['File name: ', filename(1:end-4)]);
handles.filename = Data.filename;
handles.meanImage = Data.meanImage;
handles.sizeImage = size(handles.meanImage,1);
handles.Cells = Data.Cells;
handles.trace = Data.trace;
handles.df_fixedF0 = Data.df_fixedF0;
handles.locs = Data.locs;
handles.numSpikes = Data.numSpikes;
handles.time = size(handles.trace,2);
handles.ITimeseries = Data.ITimeseries;
handles.n = 1;
handles.load = 1;
if (sum(strcmp(fieldnames(Data), 'synchrony')) == 1)
	handles.synch = Data.synchrony;
end

if (sum(strcmp(fieldnames(Data), 'cellArea')) == 1)
	handles.cellArea = Data.cellArea;
end

plotCells(handles);
plotPeaks(handles);
display('finished loading');
guidata(hObject,handles);
