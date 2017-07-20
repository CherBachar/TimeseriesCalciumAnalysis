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

% Last Modified by GUIDE v2.5 15-Jun-2017 15:33:38

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
handles.Image_fullfilepath = [filepath,filename];

set(handles.imageName, 'String', ['File name: ', filename(1:end-4)]);

I = double(imread(handles.Image_fullfilepath));
sizeI = size(I);
imgInfo = imfinfo(handles.Image_fullfilepath);
numLayers = length(imgInfo);

%if you have already analysed before or loaded data
% if numLayers > 1
%     ITimeseries = zeros(sizeI(1),sizeI(2),numLayers);
%     for time = 1:numLayers
%         ITimeseries(:,:,time) = double(imread(handles.Image_fullfilepath, time));
%     end
%     I = mean(ITimeseries,3);
%     handles.ITimeseries = ITimeseries;
% end
% 

axes(handles.axes1);
imagesc(I); colormap(gray); title('original image'); axis off;
axes(handles.axes2);
imagesc(I); colormap(gray); title('original image'); axis off;

handles.sizeImage = size(I,1);
handles.meanImage = I;
handles.load = 0;
if sum(strcmp(fieldnames(handles), 'ITimeseries')) == 1
    handles = rmfield(handles,'ITimeseries');
end
guidata(hObject,handles);


% --- Executes on button press in SegmentCells.
function SegmentCells_Callback(hObject, eventdata, handles)
% hObject    handle to SegmentCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%if you have already analysed before or loaded data
%if you have already analysed before or loaded data

%% load timeseries image
if sum(strcmp(fieldnames(handles), 'ITimeseries')) == 1
    ITimeseries =handles.ITimeseries;
    numLayers = size(ITimeseries,3); 
%if you didn't analyse or load
else
    uiwait(msgbox('Load the corresponding time-series tif stack'));
    [filename, filepath] = uigetfile('*.tif*');
    handles.imageStack_fullfilepath = [filepath,filename];
    imgInfo = imfinfo(handles.imageStack_fullfilepath);
    numLayers = length(imgInfo);
    for time = 1:numLayers
        ITimeseries(:,:,time) = double(imread(handles.imageStack_fullfilepath, time));
    end
end

handles.ITimeseries = ITimeseries;
%% segmentation of all cells
[cellLocations] = detectPoints(handles.meanImage, handles);
[Cells] = segmentCells(handles.meanImage, handles, cellLocations);
handles.Cells = Cells;


%% get cell masks
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
end

binaryImageResized = imresize(binaryImage, [sizeITs sizeITs]);
binaryNeuropilResized = imresize(binaryImageNeuropil, [sizeITs sizeITs]);

%% get active cells
%get variance through time series image
ITimeseriesSTD = std(ITimeseries,0, 3);

handles.ITimeseriesSTD = ITimeseriesSTD;
handles.maskResized = maskResized;

handles = findActiveCells(handles);



%% Plot cells
plotCells(handles);
guidata(hObject,handles);


% --- Executes on button press in Analyse.
function Analyse_Callback(hObject, eventdata, handles)
% hObject    handle to Analyse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Get cells and neuropil

%Cells
%Cells = handles.Cells;
ITimeseriesSTD =handles.ITimeseriesSTD;
ITimeseries = handles.ITimeseries;
time = size(ITimeseries,3);
activeCells = handles.activeCells;


traceNeuropil = zeros(time,1);

AllCellMask = zeros(size(ITimeseriesSTD,1),size(ITimeseriesSTD,1));
%Cells 
for l = 1:length(activeCells)
    %Get cell mask
    tempMaskCell = zeros(size(ITimeseriesSTD,1),size(ITimeseriesSTD,1));
    index = activeCells(l).PixelIdxList;
    tempMaskCell(index)=1;
    AllCellMask(index)=1;
    activeCellsMask(:,:,l) = tempMaskCell;
end
%Neuropil
neuropilMask = ITimeseriesSTD > mean(mean(ITimeseriesSTD));
neuropilMask = (neuropilMask - AllCellMask) >0;
handles.neuropilMask = neuropilMask;

%plots
figure;imagesc(AllCellMask);colormap(gray);
title('Cell mask image');
axis off
figure;imagesc(neuropilMask);colormap(gray);
title('Neuropil mask image');
axis off


%% get traces
for c = 1:length(activeCells)
    maskTempCell = logical(activeCellsMask(:,:,c)); 
    if c == 1
        maskTempNeuropil = logical(neuropilMask); 
    end
    for t = 1:time
        tempIT1 = ITimeseries(:,:,t);
        trace(c,t) = mean(mean(tempIT1(maskTempCell)));
        sumPixels(c, t) = sum(sum(tempIT1(maskTempCell)));
        if c == 1
            traceNeuropil(t) = mean(mean(tempIT1(maskTempNeuropil)));
            %sumPixels_neuropil(t) = sum(sum(tempIT1(maskTempNeuropil)));
        end
    end
    disp(['Cell: ', num2str(c)]);
    
end

handles.trace = trace;
handles.time = time;
handles.sumPixels = sumPixels;
handles.traceNeuropil = traceNeuropil;
handles.activeCellsMask = activeCellsMask;
%% Get df_F0 trace
[df_F0] = df_F0Trace(trace, handles);
handles.df_F0 = df_F0;
% [df_F0_Neuropil] = df_F0Trace(traceNeuropil,handles);
% handles.df_F0_Neuropil = df_F0_Neuropil;
%% Find synchrony 
[handles] = activityIntegral(handles);
handles.n = 1;
plotTrace(handles);

figure;imagesc(handles.corrMatrix);colormap(gray);
colorbar;
title('Correlation matrix between active cells');


%% Plot line scan
plotLineScan(handles, activeCellsMask);

disp('Finished analysis');

guidata(hObject,handles);

function CellNum_Callback(hObject, eventdata, handles)
% hObject    handle to CellNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CellNum as text
%        str2double(get(hObject,'String')) returns contents of CellNum as a double
n = str2double(get(hObject,'String'));
axes(handles.axes3);
df_F0 = handles.df_F0;
%trace = handles.trace;

%Plot in gui
plot(1:1:handles.time,df_F0(n,:));
title(['Number of active cells detected: ', num2str(length(handles.activeCells))]);
% if ~isempty(locs{n})
%     locs = handles.locs;
%     hold on 
%     plot(locs{n},df_F0(n,locs{n}), 'r*');
%     hold off
% end
handles.n = n;
plotTrace(handles);

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
    %plotCells(handles);
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
    %plotCells(handles);
end

guidata(hObject,handles);


% --- Executes on button press in addActiveCells.
function addActiveCells_Callback(hObject, eventdata, handles)
% hObject    handle to addActiveCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Read cell locations
activeCellSize = 20;
Cells = handles.activeCells;
ITimeseriesSTD = handles.ITimeseriesSTD;
imageSize = size(ITimeseriesSTD,1);
newCellsImage = zeros(imageSize,imageSize);

%% Get new cells
% [xAdd,yAdd] = ginput;
%Add a rectangle
rect = getrect;


for i = 1:size(rect,1)
%     centroid = [xAdd(i), yAdd(i)];
%     x1 = round(centroid(1) - round(activeCellSize/2));
%     x2 = round(centroid(1) + round(activeCellSize/2));
%     y1 = round(centroid(2) - round(activeCellSize/2));
%     y2 = round(centroid(2) + round(activeCellSize/2));
    x1 = rect(1); %xmin
    x2 = x1 + rect(3); %width
    y1 = rect(2);
    y2 = y1 + rect(4);
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
    
    newCellsImage(y1:y2,x1:x2) = 1;

end

newCells = regionprops(newCellsImage>0, 'Area', 'Centroid', 'Eccentricity', 'PixelIdxList', 'BoundingBox');

Cells = [Cells; newCells];

handles.activeCells = Cells;
guidata(hObject,handles);


% --- Executes on button press in removeActiveCells.
function removeActiveCells_Callback(hObject, eventdata, handles)
% hObject    handle to removeActiveCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Cells = handles.activeCells;
for c = 1:length(Cells)
    cellLocations(c,:) = Cells(c).Centroid;
end

%Get cells to remove
[xRemove,yRemove] = ginput;
if ~isempty(xRemove)
    [~, removedCells] = removeCells(cellLocations, xRemove, yRemove);
    
    Cells(removedCells) = [];
    handles.activeCells = Cells;
    %plotCells(handles);
end

guidata(hObject,handles);



% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Data.filename = handles.filename;
Data.meanImage = handles.meanImage;
Data.Cells = handles.Cells;
Data.activeCells = handles.activeCells;
Data.percentActive = handles.percentActive;

if (sum(strcmp(fieldnames(handles), 'df_F0')) == 1)
    Data.ITimeseries = handles.ITimeseries;
    Data.trace = handles.trace;
    Data.df_F0 = handles.df_F0;
    Data.activity = handles.activity;
    Data.corrMatrix = handles.corrMatrix;
    Data.ITimeseriesSTD = handles.ITimeseriesSTD;
    Data.traceNeuropil = handles.traceNeuropil;
    Data.corrPerCell = handles.corrPerCell;
    Data.corrAll = handles.corrAll;

end
save([handles.filename, '.mat'], 'Data');


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
handles.n = 1;
handles.activeCells = Data.activeCells;
handles.percentActive = Data.percentActive;


if (sum(strcmp(fieldnames(Data), 'df_F0')) == 1)
    handles.trace = Data.trace;
    handles.df_F0 = Data.df_F0;
    handles.time = size(handles.trace,2);
    handles.ITimeseries = Data.ITimeseries;
    handles.ITimeseriesSTD = Data.ITimeseriesSTD;
    handles.traceNeuropil = Data.traceNeuropil;
    handles.activity= Data.activity;
    handles.corrMatrix = Data.corrMatrix;
    plotTrace(handles);
end
if (sum(strcmp(fieldnames(Data), 'corrAll')) == 1)
    handles.corrPerCell = Data.corrPerCell;
    handles.corrAll = Data.corrAll;
end

plotCells(handles);

display('finished loading');
guidata(hObject,handles);


% --- Executes on button press in plotCells.
function plotCells_Callback(hObject, eventdata, handles)
% hObject    handle to plotCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotCells(handles);

