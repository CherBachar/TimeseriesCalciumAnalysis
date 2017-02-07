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

% Last Modified by GUIDE v2.5 07-Feb-2017 10:05:05

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
handles.meanImage_fullfilepath = [filepath,filename];
I = double(imread(handles.meanImage_fullfilepath));
axes(handles.axes1);
imagesc(I); colormap(gray); title('original image'); axis off;
axes(handles.axes2);
imagesc(I); colormap(gray); title('original image'); axis off;

handles.sizeImage = size(I,1);
handles.meanImage = I;
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
uiwait(msgbox('Load the corresponding time-series tif stack'));
[filename, filepath] = uigetfile('*.tif*');
handles.imageStack_fullfilepath = [filepath,filename];
imgInfo = imfinfo(handles.imageStack_fullfilepath);
numLayers = length(imgInfo);
Cells = handles.Cells;

for time = 1:numLayers
    ITimeseries(:,:,time) = double(imread(handles.imageStack_fullfilepath, time));
end
ITimeseries = ITimeseries;
sizeITs = size(ITimeseries,1);

binaryImage = zeros(handles.sizeImage);
for l = 1:length(Cells)
    %cellLocations(l,:) = Cells(l).Centroid;
    index = Cells(l).PixelIdxList;
    binaryImage(index)=1;
    mask(:,:,l) = binaryImage;
    maskResized(:,:,l) = imresize(mask(:,:,l), [sizeITs sizeITs]);
end


binaryImageResized = imresize(binaryImage, [sizeITs sizeITs]);
%[binaryImageLabelledResized, numCells] = bwlabel(binaryImage2, handles.con);

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
pause(1);


disp('in progress');
for c = 1:length(Cells)
    maskTemp = logical(maskResized(:,:,c)); 
    %mask3D = repmat(mask,[1,1,numLayers]); %Less efficient
    %maskOnly = mask3D .* ITimeseries;
    for time = 1:numLayers
        temp = ITimeseries(:,:,time);
        trace(c,time) = mean(mean(temp(maskTemp)));
    end    
    disp(['Cell: ', num2str(c)]);
end

handles.trace = trace;
handles.time = numLayers;
handles.ITimeseries=ITimeseries;


% Detect peaks on traces
[handles] = detectPeaks(trace, handles);
axes(handles.axes3);
df_fixedF0 = handles.df_fixedF0;
locs = handles.locs;
plot(1:1:handles.time,df_fixedF0(1,:));
title(['Number of cells detected: ', num2str(length(Cells))]);
hold on
plot(locs{1},df_fixedF0(1,locs{1}), 'r*');
hold off

disp('Finished');
guidata(hObject,handles);

% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Data.meanImage_fullfilepath = handles.meanImage_fullfilepath;
Data.imageStack_fullfilepath = handles.imageStack_fullfilepath;
Data.meanImage = handles.meanImage;
Data.Cells = handles.Cells;
Data.ITimeseries = handles.ITimeseries;
Data.trace = handles.trace;
Data.df_fixedF0 = handles.df_fixedF0;
Data.locs = handles.locs;
Data.numSpikes = handles.numSpikes;
save('Data.mat', 'Data');



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
locs = handles.locs;
plot(1:1:handles.time,df_fixedF0(n,:));
title(['Number of cells detected: ', num2str(length(locs))]);
hold on
plot(locs{n},df_fixedF0(n,locs{n}), 'r*');
hold off

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


% --- Executes on button press in RemovePeaks.
function RemovePeaks_Callback(hObject, eventdata, handles)
% hObject    handle to RemovePeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
