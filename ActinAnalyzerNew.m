function varargout = ActinAnalyzerNew(varargin)
% ACTINANALYZERNEW M-file for ActinAnalyzerNew.fig
%      ACTINANALYZERNEW, by itself, creates a new ACTINANALYZERNEW or raises the existing
%      singleton*.
%
%      H = ACTINANALYZERNEW returns the handle to a new ACTINANALYZERNEW or the handle to
%      the existing singleton*.
%
%      ACTINANALYZERNEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ACTINANALYZERNEW.M with the given input arguments.
%
%      ACTINANALYZERNEW('Property','Value',...) creates a new ACTINANALYZERNEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ActinAnalyzerNew_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ActinAnalyzerNew_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ActinAnalyzerNew

% Last Modified by GUIDE v2.5 22-Aug-2014 12:02:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ActinAnalyzerNew_OpeningFcn, ...
                   'gui_OutputFcn',  @ActinAnalyzerNew_OutputFcn, ...
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


% --- Executes just before ActinAnalyzerNew is made visible.
function ActinAnalyzerNew_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ActinAnalyzerNew (see VARARGIN)

axes(handles.axes1); axis off

% Choose default command line output for ActinAnalyzerNew
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = ActinAnalyzerNew_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in OpenImage.
function OpenImage_Callback(hObject, eventdata, handles)
% hObject    handle to OpenImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[File,Path]=uigetfile('*.tif','Choose an image');
Image=imread([Path File],'tif');
set(handles.ViewImage,'UserData',Image);
axes(handles.axes1)
imshow(Image,[]); 
axes(handles.axes1)
imshow(Image,[]); 
for i=1:4
    A=max(size(File));
    File(A)=[];
end
Identifier=[Path File];
OpenData{1}=Image;
OpenData{2}=Identifier;
set(handles.OpenImage,'UserData',OpenData);
set(handles.SelectRectangle,'UserData',Image);
handles.ThresholdText=0;

% Update handles structure
guidata(hObject, handles);


function PixelSpacing_Callback(hObject, eventdata, handles)
% hObject    handle to PixelSpacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PixelSpacing as text
%        str2double(get(hObject,'String')) returns contents of PixelSpacing as a double


% --- Executes during object creation, after setting all properties.
function PixelSpacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PixelSpacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ChooseBackgroundRegion.
function ChooseBackgroundRegion_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseBackgroundRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CurrentZ=1;
axes(handles.axes1)
Image=get(handles.ViewImage,'UserData');
PixelWidth=str2double(get(handles.PixelSpacing,'String'))./1000;
PixelArea=PixelWidth.^2;
hold on; axis off;
Background=roipoly;
      
MeanBackgroundIntensity=sum(sum(sum(double(Image).*double(Background))))./sum(sum(sum(Background)));
BackgroundValues=double(Image(Background==1));
Area=sum(sum(Background(:,:,1))).*PixelArea;
BackgroundData(1,1)=Area;
BackgroundData(2,1)=MeanBackgroundIntensity;
BackgroundData(3,1)=std(BackgroundValues);
BackgroundData(4,1)=median(BackgroundValues);
BackgroundData(5,1)=entropy(uint8(BackgroundValues));
BackgroundData(6,1)=0; %PerimeterAreaIndex;
BackgroundData(7,1)=0; %Perimeter.*PixelWidth;
BackgroundData(8,1)=0; %ShapeStats(1).MajorAxisLength;
BackgroundData(9,1)=0; %ShapeStats(1).Orientation;
BackgroundData(10,1)=0; % ShapeStats(1).Eccentricity;
Data{1}=Background;
Data{2}='N/A';
Data{3}=BackgroundData;

    Min5=ordfilt2(Image(:,:,CurrentZ),1,ones(5,5));
    Filter=fspecial('gaussian',15,5);
    MinB=imfilter(Min5,Filter);
    NewImage=Image-MinB;
    NewImage=mat2gray(NewImage);
    Xf=fspecial('gaussian',5,0.8);
    Apu=deconvlucy(NewImage,Xf);
    Helping{1}=NewImage;
    Helping{2}=Apu;
    set(handles.MeshCheckbox,'UserData',Helping);

set(handles.ChooseBackgroundRegion,'UserData',Data);
handles.BackgroundText = BackgroundData;

% Update handles structure
guidata(hObject, handles);


function VesicleDiameter_Callback(hObject, eventdata, handles)
% hObject    handle to VesicleDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function VesicleDiameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VesicleDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CentralPosPercent_Callback(hObject, eventdata, handles)
% hObject    handle to CentralPosPercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CentralPosPercent as text
%        str2double(get(hObject,'String')) returns contents of CentralPosPercent as a double


% --- Executes during object creation, after setting all properties.
function CentralPosPercent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CentralPosPercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CentralHoleNo_Callback(hObject, eventdata, handles)
% hObject    handle to CentralHoleNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CentralHoleNo as text
%        str2double(get(hObject,'String')) returns contents of CentralHoleNo as a double


% --- Executes during object creation, after setting all properties.
function CentralHoleNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CentralHoleNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CentralMeanHoleArea_Callback(hObject, eventdata, handles)
% hObject    handle to CentralMeanHoleArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CentralMeanHoleArea as text
%        str2double(get(hObject,'String')) returns contents of CentralMeanHoleArea as a double


% --- Executes during object creation, after setting all properties.
function CentralMeanHoleArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CentralMeanHoleArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CentralVesiclePenetrable_Callback(hObject, eventdata, handles)
% hObject    handle to CentralVesiclePenetrable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CentralVesiclePenetrable as text
%        str2double(get(hObject,'String')) returns contents of CentralVesiclePenetrable as a double


% --- Executes during object creation, after setting all properties.
function CentralVesiclePenetrable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CentralVesiclePenetrable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PeripheralPosPercent_Callback(hObject, eventdata, handles)
% hObject    handle to PeripheralPosPercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PeripheralPosPercent as text
%        str2double(get(hObject,'String')) returns contents of PeripheralPosPercent as a double


% --- Executes during object creation, after setting all properties.
function PeripheralPosPercent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeripheralPosPercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PeripheralHoleNo_Callback(hObject, eventdata, handles)
% hObject    handle to PeripheralHoleNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PeripheralHoleNo as text
%        str2double(get(hObject,'String')) returns contents of PeripheralHoleNo as a double


% --- Executes during object creation, after setting all properties.
function PeripheralHoleNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeripheralHoleNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PeripheralMeanHoleArea_Callback(hObject, eventdata, handles)
% hObject    handle to PeripheralMeanHoleArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PeripheralMeanHoleArea as text
%        str2double(get(hObject,'String')) returns contents of PeripheralMeanHoleArea as a double


% --- Executes during object creation, after setting all properties.
function PeripheralMeanHoleArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeripheralMeanHoleArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PeripheralVesiclePenetrable_Callback(hObject, eventdata, handles)
% hObject    handle to PeripheralVesiclePenetrable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PeripheralVesiclePenetrable as text
%        str2double(get(hObject,'String')) returns contents of PeripheralVesiclePenetrable as a double


% --- Executes during object creation, after setting all properties.
function PeripheralVesiclePenetrable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeripheralVesiclePenetrable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function ThresholdSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
Tmin=get(hObject,'Min');
Tmax=get(hObject,'Max');
Tpos=get(hObject,'Value');
ThresholdM=2^Tpos;
set(handles.ThresholdMultiplierValue,'String',mat2str(ThresholdM));

ModeT=get(handles.BlobCheckbox,'Value');
if ModeT==1
    Image=get(handles.ViewImage,'UserData');
elseif ModeT==0
    Image=get(handles.SelectRectangle,'UserData');
    Helping=get(handles.MeshCheckbox,'UserData');
    NewImage=Helping{1};
    Apu=Helping{2};
    axes(handles.axes1); hold off;
    imshow(NewImage);   
end

Centrals=handles.ChooseCentralZone;
CentralMask=Centrals{1};
PeripheralMask=handles.StoreForRingMask;

TotalMask=CentralMask+PeripheralMask;

if ModeT==1
    Mode=handles.ThresholdText;
    if Mode==0
        Mask=ones(size(Image));
    elseif Mode==1
        Mask=CentralMask;
    elseif Mode==2
        Mask=PeripheralMask;
    end
    Image=mat2gray(Image);
    Values=Image(Mask==1);
    Level=graythresh(Values);
    Threshold=ThresholdM.*Level;
    Positive=im2bw(Image,Threshold);
    axes(handles.axes1)
    imshow(Positive,[]);
    if Mode==1
        set(handles.ThresholdMultiplierValue,'UserData',Positive);
    elseif Mode==2
        set(handles.ThresholdSlider,'UserData',Positive);
    end   
elseif ModeT==0
    Values=Apu(TotalMask==1);
    Level=graythresh(Values);
    Threshold=ThresholdM.*Level;
    Positive=im2bw(Apu,Threshold);
    Positive=TotalMask.*Positive;
    LabelIm=bwlabel(Positive);
    AreaStats=regionprops(LabelIm,'Area');
    idx=find('AreaStats.Area'>1);
    Positive=ismember(LabelIm,idx);
    axes(handles.axes1)
    imshow(Positive,[]);
    Pos1=Positive.*CentralMask;
    Pos2=Positive.*PeripheralMask;
    set(handles.ThresholdMultiplierValue,'UserData',Pos1);
    set(handles.ThresholdSlider,'UserData',Pos2);
end

% Update handles structure
guidata(hObject, handles);


function ThresholdSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Min',-4,'Max',3)
set(hObject,'Value',0)

% Update handles structure
guidata(hObject, handles);


function ThresholdMultiplierValue_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdMultiplierValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ThresholdMultiplierValue as text
%        str2double(get(hObject,'String')) returns contents of ThresholdMultiplierValue as a double
ThresholdM=str2double(get(hObject,'String'));
SliderPos=log2(ThresholdM);
set(handles.ThresholdSlider,'Value',SliderPos);
ModeT=get(handles.BlobCheckbox,'Value');
if ModeT==1
    Image=get(handles.ViewImage,'UserData');
elseif ModeT==0
    Image=get(handles.SelectRectangle,'UserData');
    Helping=get(handles.MeshCheckbox,'UserData');
    NewImage=Helping{1};
    Apu=Helping{2};
    axes(handles.axes1); hold off;
    imshow(NewImage);   
end
Centrals=handles.ChooseCentralZone;
CentralMask=Centrals{1};
PeripheralMask=handles.StoreForRingMask;
TotalMask=CentralMask+PeripheralMask;

if ModeT==1
    Mode=handles.ThresholdText;
    if Mode==0
        Mask=ones(size(Image));
    elseif Mode==1
        Mask=CentralMask;
    elseif Mode==2
        Mask=PeripheralMask;
    end
    Image=mat2gray(Image);
    Values=Image(Mask==1);
    Level=graythresh(Values);
    Threshold=ThresholdM.*Level;
    Positive=im2bw(Image,Threshold);
    axes(handles.axes1)
    imshow(Positive,[]);
    if Mode==1
        set(handles.ThresholdMultiplierValue,'UserData',Positive);
    elseif Mode==2
        set(handles.ThresholdSlider,'UserData',Positive);
    end   
elseif ModeT==0
    Values=Apu(TotalMask==1);
    Level=graythresh(Values);
    Threshold=ThresholdM.*Level;
    Positive=im2bw(Apu,Threshold);
    Positive=TotalMask.*Positive;
    LabelIm=bwlabel(Positive);
    AreaStats=regionprops(LabelIm,'Area');
    idx=find('AreaStats.Area'>1);
    Positive=ismember(LabelIm,idx);
    axes(handles.axes1)
    imshow(Positive,[]);
    Pos1=Positive.*CentralMask;
    Pos2=Positive.*PeripheralMask;
    set(handles.ThresholdMultiplierValue,'UserData',Pos1);
    set(handles.ThresholdSlider,'UserData',Pos2);
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ThresholdMultiplierValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThresholdMultiplierValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AnalyzeMesh.
function AnalyzeMesh_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Image=get(handles.ViewImage,'UserData');

ImageSize=size(Image);
VesicleR=0.5*str2double(get(handles.VesicleDiameter,'String'))./1000;
PixelWidth=str2double(get(handles.PixelSpacing,'String'))./1000;
VesicleR=VesicleR./PixelWidth;
PixelArea=PixelWidth.^2;
Centrals=handles.ChooseCentralZone;
Peripherals=handles.ChooseCellPerimeter;
CentralMask=Centrals{1};
RingMask=handles.StoreForRingMask;
CentralPositive=get(handles.ThresholdMultiplierValue,'UserData');
RingPositive=get(handles.ThresholdSlider,'UserData');

TotalMask=CentralMask+RingMask; %NEW LINE 25th May 2011
TotalMask=imfill(TotalMask,'holes');
TotalMaskLabel=bwlabel(TotalMask,8); %NEW LINE 25th May 2011
TotalMaskStats=regionprops(TotalMaskLabel,'Area','Centroid'); %NEW LINE 25th May 2011
LTMS=length([TotalMaskStats.Area]); %NEW LINE 25th May 2011
if LTMS==1 %NEW LINE 25th May 2011
    CentroidXY=TotalMaskStats(1).Centroid; %NEW LINE 25th May 2011
else %NEW LINE 25th May 2011
    for q=1:LTMS %NEW LINE 25th May 2011
        Cent(q,:)=TotalMaskStats(q).Centroid; %NEW LINE 25th May 2011
        Are(q,:)=TotalMaskStats(q).Area; %NEW LINE 25th May 2011
    end %NEW LINE 25th May 2011
    [~,MI]=max(Are); %NEW LINE 25th May 2011
    CentroidXY=Cent(MI,:); %NEW LINE 25th May 2011
    
end %NEW LINE 25th May 2011

Data=zeros(5,2);
Data(1,1)=sum(sum(CentralPositive))./(sum(sum(CentralMask)));
set(handles.CentralPosPercent,'String',mat2str(round(100*Data(1,1)),-1));
Data(1,2)=sum(sum(RingPositive))./(sum(sum(RingMask)));
set(handles.PeripheralPosPercent,'String',mat2str(round(100*Data(1,2)),-1));

CentralDistanceMask=(1-CentralMask)+CentralPositive;
RingDistanceMask=(1-RingMask)+RingPositive;
CentralTransform=bwdist(CentralDistanceMask);
CentralPenet=zeros(ImageSize);
CentralPenet(CentralTransform>VesicleR)=1;
CutOffDistanceTransforms{1}=CentralTransform.*CentralPenet.*PixelWidth.*1000;
Data(5,1)=sum(sum(CentralPenet))./(sum(sum(CentralMask)));
set(handles.CentralVesiclePenetrable,'String',mat2str((100.*Data(5,1)),-1));

RingTransform=bwdist(RingDistanceMask);
RingPenet=zeros(ImageSize);
RingPenet(RingTransform>VesicleR)=1;
CutOffDistanceTransforms{2}=RingTransform.*RingPenet.*PixelWidth.*1000;
Data(5,2)=sum(sum(RingPenet))./(sum(sum(RingMask)));
set(handles.PeripheralVesiclePenetrable,'String',mat2str(round(100*Data(5,2)),-1));

CutOffDistanceTransforms{3}=(RingTransform.*RingPenet+CentralTransform.*CentralPenet).*PixelWidth.*1000;

MaxVes=VesicleR.*1000.*PixelWidth+255;
VesicleRnm=VesicleR.*1000.*PixelWidth;
CutOffDistanceTransforms{3}=0.5*((CutOffDistanceTransforms{3}-VesicleRnm)+abs(CutOffDistanceTransforms{3}-VesicleRnm));
CutOffDistanceTransforms{3}=CutOffDistanceTransforms{3}./MaxVes;
max(max(CutOffDistanceTransforms{3}));
CutOffDistanceTransforms{3}(CutOffDistanceTransforms{3}>1)=1;
CutOffDistanceTransforms{3}=round(255.*CutOffDistanceTransforms{3});

set(handles.ViewPenetration,'UserData',CutOffDistanceTransforms);

CentralHoles=CentralMask.*(1-CentralPositive);
RingHoles=RingMask.*(1-RingPositive);

LinInd=find(CentralPenet==1);  %NEW LINE 25th May 2011
[PenetY,PenetX]=ind2sub(size(CentralPenet),LinInd);  %NEW LINE 25th May 2011
PenetY=double(PenetY); %NEW LINE 25th May 2011
PenetX=double(PenetX); %NEW LINE 25th May 2011
LinInd=find(CentralMask==1);  %NEW LINE 25th May 2011
[CentralY,CentralX]=ind2sub(size(CentralMask),LinInd);  %NEW LINE 25th May 2011
CentralY=double(CentralY); %NEW LINE 25th May 2011
CentralX=double(CentralX); %NEW LINE 25th May 2011
LinInd=find(CentralHoles==1);  %NEW LINE 25th May 2011
[HolesY,HolesX]=ind2sub(size(CentralHoles),LinInd);  %NEW LINE 25th May 2011
HolesY=double(HolesY); %NEW LINE 25th May 2011
HolesX=double(HolesX); %NEW LINE 25th May 2011
PenetDist=PixelWidth.*sqrt((PenetX-CentroidXY(1,1)).^2+(PenetY-CentroidXY(1,2)).^2); %NEW LINE 25th May 2011
CentralDist=PixelWidth.*sqrt((CentralX-CentroidXY(1,1)).^2+(CentralY-CentroidXY(1,2)).^2); %NEW LINE 25th May 2011
HolesDist=PixelWidth.*sqrt((HolesX-CentroidXY(1,1)).^2+(HolesY-CentroidXY(1,2)).^2); %NEW LINE 25th May 2011
DistData(1,1)=mean(PenetDist); %NEW LINE 25th May 2011
DistData(1,2)=std(PenetDist); %NEW LINE 25th May 2011
DistData(1,3)=min(PenetDist); %NEW LINE 25th May 2011
DistData(2,1)=mean(CentralDist); %NEW LINE 25th May 2011
DistData(2,2)=std(CentralDist); %NEW LINE 25th May 2011
DistData(2,3)=min(CentralDist); %NEW LINE 25th May 2011
DistData(3,1)=mean(HolesDist); %NEW LINE 25th May 2011
DistData(3,2)=std(HolesDist); %NEW LINE 25th May 2011
DistData(3,3)=min(HolesDist); %NEW LINE 25th May 2011

mask_panel = figure;
set(mask_panel, 'Position', [50, 50, 750, 250]);
subaxis(1,3,1,'Padding',0,'Spacing',0,'Margin',0); imshow(TotalMask);
subaxis(1,3,2,'Padding',0,'Spacing',0,'Margin',0); imshow(CentralHoles);
subaxis(1,3,3,'Padding',0,'Spacing',0,'Margin',0); imshow(RingHoles);
pause(2);
close(mask_panel)

HoleIm{1}=CentralHoles+RingHoles;
HoleIm{2}=[];

CentralLabel=bwlabel(CentralHoles,8);
RingLabel=bwlabel(RingHoles,8);

HoleStats=regionprops(CentralLabel,'Area');
k=1;
MaxI=max(size(HoleStats));
Areas=[];
idx=find([HoleStats.Area]>2);
HoleIm{3}=ismember(CentralLabel,idx);

for i=1:MaxI
    if HoleStats(i).Area>2
        Areas(k,1)=PixelArea.*HoleStats(i).Area;
        k=k+1;
    end
end
CentralHoleAreas=Areas;
Data(2,1)=max(size(Areas));
set(handles.CentralHoleNo,'String',mat2str(Data(2,1)));
Data(3,1)=mean(Areas);
set(handles.CentralMeanHoleArea,'String',mat2str(Data(3,1)));
Data(4,1)=median(Areas);
set(handles.CentralMedianHoleArea,'String',mat2str(Data(4,1)));
Data(6,1)=std(Areas);
Data(7,1)=Data(6,1)./sqrt(Data(2,1));
Data(8,1)=prctile(Areas,5);
Data(9,1)=prctile(Areas,95);
Data(10,1)=prctile(Areas,25);
Data(11,1)=prctile(Areas,75);

HoleStats=regionprops(RingLabel,'Area','Extrema','Extent');
k=1;
MaxI=max(size(HoleStats));
Areas=[];
idx=find([HoleStats.Area]>2);
HoleIm{2}=ismember(RingLabel,idx);

Ux=zeros(8,8);
Uy=zeros(8,8);
Vx=zeros(8,8);
Vy=zeros(8,8);
for i=1:MaxI
    Mat1=HoleStats(i).Extrema;
    Mat1X=Mat1(:,1);
    Mat1Y=Mat1(:,2);
    Mat2X=transpose(Mat1X);
    Mat2Y=transpose(Mat1Y);
    for z=1:8
       Ux(:,z)=Mat1X;
       Uy(:,z)=Mat1Y;
       Vx(z,:)=Mat2X;
       Vy(z,:)=Mat2Y;
    end
    Distance=sqrt((Ux-Vx).^2+(Uy-Vy).^2);
    MaxD=max(max(Distance));
    TestValue=HoleStats(i).Area;
    if (TestValue>0.1)&&(HoleStats(i).Extent>0.1)
    if HoleStats(i).Area>2
        Areas(k,1)=PixelArea.*HoleStats(i).Area;
        idx2(k)=i;
        k=k+1;
    end
    end
end
HoleIm{4}=ismember(RingLabel,idx2);
HoleIm{2}=HoleIm{2}-HoleIm{4};

handles.ViewRingHoles=HoleIm;
PeripheralHoleAreas=Areas;
Data(2,2)=max(size(Areas));
set(handles.PeripheralHoleNo,'String',mat2str(Data(2,2)));
Data(3,2)=mean(Areas);
set(handles.PeripheralMeanHoleArea,'String',mat2str(Data(3,2)));
Data(4,2)=median(Areas);
set(handles.PeripheralMedianHoleArea,'String',mat2str(Data(4,2)));
Data(6,2)=std(Areas);
Data(7,2)=Data(6,2)./sqrt(Data(2,2));
Data(8,2)=prctile(Areas,5);
Data(9,2)=prctile(Areas,95);
Data(10,2)=prctile(Areas,25);
Data(11,2)=prctile(Areas,75);

OpenData=get(handles.OpenImage,'UserData');
File=OpenData{2};
File=[File '_central_thresholded.tif'];
imwrite(CentralPositive,File,'tif');
File=OpenData{2};
File=[File '_peripheral_holes.txt'];
dlmwrite(File,PeripheralHoleAreas,'\t');
File=OpenData{2};
File=[File '_central_holes.txt'];
dlmwrite(File,CentralHoleAreas,'\t');
File=OpenData{2};
File=[File '_peripheral_thresholded.tif'];
imwrite(RingPositive,File,'tif');
Data(5,:)=Data(5,:).*100;
Data(1,:)=Data(1,:).*100;
set(handles.AnalyzeMesh,'UserData',Data);

HoleSizePC{1}=PeripheralHoleAreas;
HoleSizePC{2}=CentralHoleAreas;
set(handles.RegionSelectionPanel,'UserData',HoleSizePC);
File=OpenData{2};
File=[File '_AnalyzeMesh_Borders.mat'];
save(File,'CutOffDistanceTransforms','CentralPositive','RingPositive','RingMask','Centrals','Peripherals');
File=OpenData{2}; %NEW LINE 25th May 2011
File=[File '_DistDataFile.mat']; %NEW LINE 25th May 2011
save(File,'CentralX','CentralY','CentroidXY','HolesX','HolesY','PenetX','PenetY','DistData','PenetDist','CentralDist','HolesDist'); %NEW LINE 25th May 2011
clear LinInd PenetY PenetX HolesX HolesY CentralX CentralY PenetDist CentralDist HolesDist %NEW LINE 25th May 2011

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function ActinAnalyzeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ActinAnalyzeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SaveMeshData.
function SaveMeshData_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMeshData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Data=get(handles.AnalyzeMesh,'UserData');
DataCell=num2cell(Data);
Top={'Central' 'Peripheral'};
DataCell=[Top; DataCell];
Edge={' '; 'Actin mesh area%'; '#holes'; 'mean hole area'; 'median hole area'; '% vesicle penetrable area'; 'hole area std'; 'hole area sem'; 'hole area 5th percentile'; 'hole area 95th percentile'; 'hole area 25th percentile'; 'hole area 75th percentile'};
DataCell=[Edge DataCell];
OpenData=get(handles.OpenImage,'UserData');
File=OpenData{2};
ModeT=get(handles.BlobCheckbox,'UserData');
if ModeT==1
    File=[File '_blob_mode_mesh_data.xls'];
    FileD=[File '_blob_mode_dist_data.xls'];
elseif ModeT==0
    File=[File '_mesh_mode_mesh_data.xls'];
    FileD=[File '_mesh_mode_dist_data.xls'];
elseif ModeT==2
    File=[File '_automatic_mode_mesh_data.xls'];
    FileD=[File '_automatic_mode_dist_data.xls'];
end
xlswrite2(File,DataCell);
HolesizePC=get(handles.RegionSelectionPanel,'UserData');
PeripheralHoleAreas=HolesizePC{1};
CentralHoleAreas=HolesizePC{2};
CellName=get(handles.CellType,'String');

File=[tempdir CellName '.mat'];
try
    %disp('ping')
    load(File);
    HoleSizeP=[HoleSizeP; PeripheralHoleAreas];
    HoleSizeC=[HoleSizeC; CentralHoleAreas];
catch
    %disp('pong')
    HoleSizeP=PeripheralHoleAreas;
    HoleSizeC=CentralHoleAreas;
end
save(File,'HoleSizeP','HoleSizeC')

File2=OpenData{2}; %NEW LINE 25th MAY 2011
File2=[File2 '_DistDataFile.mat']; %NEW LINE 25th MAY 2011
load(File2) %NEW LINE 25th MAY 2011
DataCell=num2cell(DistData); %NEW LINE 25th MAY 2011 
Top={'mean' 'std' 'min'}; %NEW LINE 25th MAY 2011 
DataCell=[Top; DataCell]; %NEW LINE 25th MAY 2011
Edge={' '; 'Centroid-Granule penetrable area distance'; 'Centroid-Central pixels distance'; 'Centroid-Holes distance'}; %NEW LINE 25th MAY 2011  % 'median hole area'; '% vesicle penetrable area'; 'hole area std'; 'hole area sem'; 'hole area 5th percentile'; 'hole area 95th percentile'; 'hole area 25th percentile'; 'hole area 75th percentile'};
DataCell=[Edge DataCell]; %NEW LINE 25th MAY 2011
File3=OpenData{2};
File3=[File3 '_DistDataFile.xls'];
xlswrite2(File3,DataCell) %NEW LINE 25th MAY 2011

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in ViewImage.
function ViewImage_Callback(hObject, eventdata, handles)
% hObject    handle to ViewImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
OpenData=get(handles.OpenImage,'UserData');
Image=OpenData{1};
axes(handles.axes1)
hold off;
imshow(Image,[]);
set(handles.ViewImage,'UserData',Image);
handles.ThresholdText=0;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in SelectRectangle.
function SelectRectangle_Callback(hObject, eventdata, handles)
% hObject    handle to SelectRectangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Image=uint8(get(handles.ViewImage,'UserData'));
size(Image)
max(max(Image))
OpenData=get(handles.OpenImage,'UserData');

axes(handles.axes1)
cla(handles.axes1)
imshow(Image,[]);
hold off;
SelectionHandle=imrect;
Position=getPosition(SelectionHandle);
Xlower=round(Position(1));
Xhigher=round(Position(1)+Position(3));
Ylower=round(Position(2));
Yhigher=round(Position(2)+Position(4));
SelectionImage=Image(Ylower:Yhigher,Xlower:Xhigher);
set(handles.ViewImage,'UserData',SelectionImage);
set(handles.SelectRectangle,'UserData',SelectionImage);
imshow(SelectionImage,[]);

    Im=get(handles.OpenMultiCh,'UserData');
    size(Im{1})
     size(Im{2})
    size(Im{3})
    for i=1:3
        Im{i}=Im{i}(Ylower:Yhigher,Xlower:Xhigher);
    end
    set(handles.OpenMultiCh,'UserData');


File=OpenData{2};
File=[File '_RectangleStack.mat'];
save(File,'Channels');

% Update handles structure
guidata(hObject, handles);


function CentralMedianHoleArea_Callback(hObject, eventdata, handles)
% hObject    handle to CentralMedianHoleArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CentralMedianHoleArea as text
%        str2double(get(hObject,'String')) returns contents of CentralMedianHoleArea as a double


% --- Executes during object creation, after setting all properties.
function CentralMedianHoleArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CentralMedianHoleArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PeripheralMedianHoleArea_Callback(hObject, eventdata, handles)
% hObject    handle to PeripheralMedianHoleArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PeripheralMedianHoleArea as text
%        str2double(get(hObject,'String')) returns contents of PeripheralMedianHoleArea as a double


% --- Executes during object creation, after setting all properties.
function PeripheralMedianHoleArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeripheralMedianHoleArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CellMeanInt_Callback(hObject, eventdata, handles)
% hObject    handle to CellMeanInt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function CellMeanInt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellMeanInt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function CellMeanInt_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to CellMeanInt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in FindPerimeters.
function FindPerimeters_Callback(hObject, eventdata, handles)
% hObject    handle to FindPerimeters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Image=get(handles.ViewImage,'UserData');

BackG=handles.BackgroundText;
Image=0.5*(single(Image)-BackG(2)+abs(single(Image)-BackG(2)));
Image=mat2gray(Image);
initialLSF=zeros(size(Image));
initialLSF(Image>graythresh(Image))=1;
initialLSF=initialLSF-0.5;
Image=round(255*single(Image));
ImageB=single(SqrtWeighed(Image));
u=RunLBFfastSingleImage(ImageB,initialLSF,19,120,handles.axes1);

    axes(handles.axes1);
    imshow(Image(:,:,1),[]);
    ContourMatrix=contour(u(:,:,1),[0 0]);
    LargeLimit=250;
    [~,LargeContours,SmallContours]=ContourOrganizerWithout(ContourMatrix,LargeLimit);
    [OuterContours,InnerContours]=PickOuterContoursSingle(LargeContours);
    axes(handles.axes1);
    hold on; axis off;
    SizeOuter=max(size(OuterContours));
    MaxOuter=[];
    MaxInner=[];
    for i=1:SizeOuter
        if max(size(MaxOuter))<max(size(OuterContours{i}));
            MaxOuter=OuterContours{i};
        end
    end
    plot(MaxOuter(:,1),MaxOuter(:,2),'m');
    SizeInner=max(size(InnerContours));
    for i=1:SizeInner
        if max(size(MaxInner))<max(size(InnerContours{i}));
            MaxInner=InnerContours{i};
        end
        
    end
    if max(size(MaxInner))>250;
        plot(MaxInner(:,1),MaxInner(:,2),'c')
    end
        DataX{1}=u(:,:,1);
        DataX{2}=MaxOuter;
        DataX{3}=MaxInner;
        DataX{4}=Image(:,:,1);

        set(handles.FindPerimeters,'UserData',DataX);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in CutU.
function CutU_Callback(hObject, eventdata, handles)
% hObject    handle to CutU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;
FindData=get(handles.FindPerimeters,'UserData');

Image=FindData{4};
imshow(Image,[]); hold on; axis off;
    u=FindData{1};
    if (max(size(FindData{2}))>250)&&(max(size(FindData{3})>250));
        plot(FindData{2}(:,1),FindData{2}(:,2),'m',FindData{3}(:,1),FindData{3}(:,2),'c');
    else
        plot(FindData{2}(:,1),FindData{2}(:,2),'m');
    end
    
FreeHandle=imfreehand;
FreehandPos=getPosition(FreeHandle);
SizeFHP=size(FreehandPos);

PosMatrix=FreehandPos(1,:);
for i=1:(SizeFHP(1)-1);
    Span=max(abs(FreehandPos(i+1,1)-FreehandPos(i,1)),abs(FreehandPos(i+1,2)-FreehandPos(i,2)));
    if Span>1
        if abs(FreehandPos(i+1,1)-FreehandPos(i,1))>0;
            InterpX=transpose(min(FreehandPos(i,1),FreehandPos(i+1,1)):(1/(Span*2)):max(FreehandPos(i,1),FreehandPos(i+1,1)));
            InterpY=interp1([FreehandPos(i,1); FreehandPos(i+1,1)], [FreehandPos(i,2); FreehandPos(i+1,2)],InterpX,'linear');
        else
            InterpY=transpose(min(FreehandPos(i,2),FreehandPos(i+1,2)):(1/(Span*2)):max(FreehandPos(i,2),FreehandPos(i+1,2)));
            InterpX=interp1([FreehandPos(i,2); FreehandPos(i+1,2)], [FreehandPos(i,1); FreehandPos(i+1,1)],InterpY,'linear');
        end
        Cop=[InterpX InterpY];
        
    else
        Cop=FreehandPos(i+1,:);
    end
    PosMatrix=[PosMatrix ;Cop];
end
FreehandPos=PosMatrix;
Coord1=[floor(FreehandPos(:,1)),floor(FreehandPos(:,2))];
Coord2=[floor(FreehandPos(:,1)),ceil(FreehandPos(:,2))];
Coord3=[ceil(FreehandPos(:,1)),floor(FreehandPos(:,2))];
Coord4=[ceil(FreehandPos(:,1)),ceil(FreehandPos(:,2))];
Coord=[Coord1; Coord2; Coord3; Coord4];

ImageSize=size(Image);
LinearIndices=sub2ind(ImageSize,Coord(:,2),Coord(:,1));
u(LinearIndices)=-100;
ContourMatrix=contour(u,[0 0]);
LargeLimit=250;
[~,LargeContours,SmallContours]=ContourOrganizerWithout(ContourMatrix,LargeLimit);
[OuterContours,InnerContours]=PickOuterContoursSingle(LargeContours);
SizeOuter=max(size(OuterContours));
MaxOuter=[];
MaxInner=[];
axes(handles.axes1);
hold off;

imshow(Image,[]); hold on; axis off;
    
for i=1:SizeOuter
    if max(size(MaxOuter))<max(size(OuterContours{i}));
        MaxOuter=OuterContours{i};
    end
end
plot(MaxOuter(:,1),MaxOuter(:,2),'m');
SizeInner=max(size(InnerContours));
for i=1:SizeInner
    if max(size(MaxInner))<max(size(InnerContours{i}));
        MaxInner=InnerContours{i};
    end
end

if max(size(MaxInner))>250;
    plot(MaxInner(:,1),MaxInner(:,2),'c');
end

    DataX{1}=u;
    DataX{2}=MaxOuter;
    DataX{3}=MaxInner;
    DataX{4}=Image;
    set(handles.FindPerimeters,'UserData',DataX);

    hold off

    % Update handles structure
guidata(hObject, handles);

% --- Executes on button press in JoinU.
function JoinU_Callback(hObject, eventdata, handles)
% hObject    handle to JoinU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1)
FindData=get(handles.FindPerimeters,'UserData');
Image=FindData{4};
imshow(Image,[]); hold on; axis off;
    u=FindData{1};
    if (max(size(FindData{2}))>250)&&(max(size(FindData{3})>250))
        plot(FindData{2}(:,1),FindData{2}(:,2),'m',FindData{3}(:,1),FindData{3}(:,2),'c')
    else
        plot(FindData{2}(:,1),FindData{2}(:,2),'m')
    end
FreeHandle=imfreehand;
FreehandPos=getPosition(FreeHandle);
SizeFHP=size(FreehandPos);

PosMatrix=FreehandPos(1,:);
for i=1:(SizeFHP(1)-1)
    Span=max(abs(FreehandPos(i+1,1)-FreehandPos(i,1)),abs(FreehandPos(i+1,2)-FreehandPos(i,2)));
    if Span>1
        if abs(FreehandPos(i+1,1)-FreehandPos(i,1))>0
            InterpX=transpose(min(FreehandPos(i,1),FreehandPos(i+1,1)):(1/(Span*2)):max(FreehandPos(i,1),FreehandPos(i+1,1)));
            InterpY=interp1([FreehandPos(i,1); FreehandPos(i+1,1)], [FreehandPos(i,2); FreehandPos(i+1,2)],InterpX,'linear');
        else
            InterpY=transpose(min(FreehandPos(i,2),FreehandPos(i+1,2)):(1/(Span*2)):max(FreehandPos(i,2),FreehandPos(i+1,2)));
            InterpX=interp1([FreehandPos(i,2); FreehandPos(i+1,2)], [FreehandPos(i,1); FreehandPos(i+1,1)],InterpY,'linear');
        end
        Cop=[InterpX InterpY];
        
    else
        Cop=FreehandPos(i+1,:);
    end
    PosMatrix=[PosMatrix ;Cop];
end
ImageSize=size(Image);
FreehandPos=PosMatrix;
Coord1=[floor(FreehandPos(:,1)),floor(FreehandPos(:,2))];
Coord2=[floor(FreehandPos(:,1)),ceil(FreehandPos(:,2))];
Coord3=[ceil(FreehandPos(:,1)),floor(FreehandPos(:,2))];
Coord4=[ceil(FreehandPos(:,1)),ceil(FreehandPos(:,2))];
Coord=[Coord1; Coord2; Coord3; Coord4];
Help=zeros(ImageSize);
LinInd=sub2ind(ImageSize,Coord(:,2),Coord(:,1));
Help(LinInd)=1;
Help=imfill(Help,4,'holes');

LinearIndices=sub2ind(ImageSize,Coord(:,2),Coord(:,1));
u(LinearIndices)=100;
u(Help==1)=100;
FindData{1}=u;
set(handles.FindPerimeters,'UserData',FindData);
ContourMatrix=contour(u,[0 0]);
LargeLimit=250;
[~,LargeContours,SmallContours]=ContourOrganizerWithout(ContourMatrix,LargeLimit);
[OuterContours,InnerContours]=PickOuterContoursSingle(LargeContours);

SizeOuter=max(size(OuterContours));
MaxOuter=[];
MaxInner=[];
axes(handles.axes1)
hold off;
    
for i=1:SizeOuter
    if max(size(MaxOuter))<max(size(OuterContours{i}));
        MaxOuter=OuterContours{i};
    end
end
plot(MaxOuter(:,1),MaxOuter(:,2),'m');
SizeInner=max(size(InnerContours));
for i=1:SizeInner
    if max(size(MaxInner))<max(size(InnerContours{i}));
        MaxInner=InnerContours{i};
    end
    
end
if max(size(MaxInner))>250
    plot(MaxInner(:,1),MaxInner(:,2),'c')
end
    DataX{1}=u;
    DataX{2}=MaxOuter;
    DataX{3}=MaxInner;
    DataX{4}=Image;
    set(handles.FindPerimeters,'UserData',DataX)

save([tempdir,'Tests.mat'],'u','OuterContours','InnerContours','MaxOuter','MaxInner');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in UseFound.
function UseFound_Callback(hObject, eventdata, handles)
% hObject    handle to UseFound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    FindData=get(handles.FindPerimeters,'UserData');
    Outer=FindData{2};
    Inner=FindData{3};
    ImageSize=size(FindData{4});
    All=[];
    Xsize=ImageSize(2);
    Ysize=ImageSize(1);
    for i=1:Xsize
        X=ones(Ysize,1).*i;
        Y=(1:1:Ysize)';
        All=[All; X Y];
    end
    CellPerim=zeros(ImageSize);
    InnerPerim=zeros(ImageSize);
    Coord1=[floor(Outer(:,1)),floor(Outer(:,2))];
    Coord2=[floor(Outer(:,1)),ceil(Outer(:,2))];
    Coord3=[ceil(Outer(:,1)),floor(Outer(:,2))];
    Coord4=[ceil(Outer(:,1)),ceil(Outer(:,2))];
    Coord=[Coord1; Coord2; Coord3; Coord4];
    LinearIndices=sub2ind(ImageSize,Coord(:,2),Coord(:,1));
    CellPerim(LinearIndices)=1;
    CellPerim=imfill(CellPerim,'holes');
    
    Coord1=[floor(Inner(:,1)),floor(Inner(:,2))];
    Coord2=[floor(Inner(:,1)),ceil(Inner(:,2))];
    Coord3=[ceil(Inner(:,1)),floor(Inner(:,2))];
    Coord4=[ceil(Inner(:,1)),ceil(Inner(:,2))];
    Coord=[Coord1; Coord2; Coord3; Coord4];
    LinearIndices=sub2ind(ImageSize,Coord(:,2),Coord(:,1));
    InnerPerim(LinearIndices)=1;
    InnerPerim=imfill(InnerPerim,'holes');
    
    PixelWidth=str2double(get(handles.PixelSpacing,'String'))./1000;
    PixelArea=PixelWidth.^2;
    
    Image=get(handles.ViewImage,'UserData');
    
    CellP{1}=CellPerim;
    BackgroundMask=1-CellPerim;
    set(handles.ViewBackground,'UserData',BackgroundMask);
    
    CellP{2}=Outer;
    MeanCellIntensity=sum(sum(double(Image).*double(CellPerim)))./sum(sum(CellPerim));
    CellValues=double(Image(CellPerim==1));
    [ConvexAreaToArea,PerimeterAreaIndex,ContourArea,ConvexArea,Perimeter]=CheckContour(Outer);
    ShapeStats=regionprops(CellPerim,'Eccentricity','MajorAxisLength','Orientation');
    Area=sum(sum(CellPerim)).*PixelArea;
    Data(1,1)=Area;
    Data(2,1)=MeanCellIntensity;
    Data(3,1)=std(CellValues);
    Data(4,1)=median(CellValues);
    Data(5,1)=entropy(uint8(CellValues));
    Data(6,1)=PerimeterAreaIndex;
    Data(7,1)=Perimeter.*PixelWidth;
    Data(8,1)=ShapeStats(1).MajorAxisLength;
    Data(9,1)=ShapeStats(1).Orientation;
    Data(10,1)=ShapeStats(1).Eccentricity;
    CellP{3}=Data;
    handles.ChooseCellPerimeter=CellP;
    
    IntP{1}=InnerPerim;
    IntP{2}=Inner;
    MeanCentralIntensity=sum(sum(double(Image).*double(InnerPerim)))./sum(sum(InnerPerim));
    CentralValues=double(Image(InnerPerim==1));
    [ConvexAreaToArea,PerimeterAreaIndex,ContourArea,ConvexArea,Perimeter]=CheckContour(Inner);
    ShapeStats=regionprops(InnerPerim,'Eccentricity','MajorAxisLength','Orientation');
    Area=sum(sum(InnerPerim)).*PixelArea;
    Data(1,1)=Area;
    Data(2,1)=MeanCentralIntensity;
    Data(3,1)=std(CentralValues);
    Data(4,1)=median(CentralValues);
    Data(5,1)=entropy(uint8(CentralValues));
    Data(6,1)=PerimeterAreaIndex;
    Data(7,1)=Perimeter.*PixelWidth;
    Data(8,1)=ShapeStats(1).MajorAxisLength;
    Data(9,1)=ShapeStats(1).Orientation;
    Data(10,1)=ShapeStats(1).Eccentricity;
    IntP{3}=Data;
    handles.ChooseCentralZone=IntP;
    RingMask=CellPerim.*(1-InnerPerim);
    handles.StoreForRingMask=RingMask; 

% Update handles structure
guidata(hObject, handles);

function CellType_Callback(hObject, eventdata, handles)
% hObject    handle to CellType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function CellType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in HoleSizeHeatmap.
function HoleSizeHeatmap_Callback(hObject, eventdata, handles)
% hObject    handle to HoleSizeHeatmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Image=get(handles.SelectRectangle,'UserData');
ImageSize=size(Image);
PixelWidth=str2double(get(handles.PixelSpacing,'String'))./1000;
PixelArea=PixelWidth.^2;
%if largest hole size is 3 sq. micrometers
LinMaxHole=3./PixelArea; %hole size in pixels

Mode=get(handles.SquareRootMode,'UserData');
if isempty(Mode)==1
    Mode=1;
end
Cmap1=colormap(gray(256));

HoleIm=handles.ViewRingHoles;

CentralHoles=HoleIm{3};
RingHoles=HoleIm{4};

CentralLabel=bwlabel(CentralHoles,8);
RingLabel=bwlabel(RingHoles,8);

HoleStats=regionprops(CentralLabel,'Area');
MaxI=max(size(HoleStats));
AreasC=[];
for i=1:MaxI
        AreasC(i,1)=HoleStats(i).Area;
        AreasC(i,2)=i;
end

HoleStats=regionprops(RingLabel,'Area');
MaxI=max(size(HoleStats));
AreasP=[];
for i=1:MaxI
        AreasP(i,1)=HoleStats(i).Area;
        AreasP(i,2)=i;
end

MaxC=max(AreasC,[],1);
MaxC=MaxC(1);
MaxP=max(AreasP,[],1);
MaxP=MaxP(1);
MaxH=max(MaxC,MaxP);
if Mode==1
    RSqMax=round(sqrt(LinMaxHole));
elseif Mode==0
    RSqMax=round(LinMaxHole);
elseif Mode==2
    RSqMax=round(log(LinMaxHole));
end
    
Cmap2=colormap(jet(255));
IndexC=zeros(size(AreasC));
IndexP=zeros(size(AreasP));
if Mode==1
    IndexC(:,1)=round(255*sqrt(AreasC(:,1))./RSqMax);
    IndexP(:,1)=round(255*sqrt(AreasP(:,1))./RSqMax);
elseif Mode==0
    IndexC(:,1)=round(255*(AreasC(:,1))./RSqMax);
    IndexP(:,1)=round(255*(AreasP(:,1))./RSqMax);
elseif Mode==2
    IndexC(:,1)=round(255*log(AreasC(:,1))./RSqMax);
    IndexP(:,1)=round(255*log(AreasP(:,1))./RSqMax);
end
IndexC(IndexC>255)=255;
IndexP(IndexP>255)=255;
IndexC=IndexC+(1-ceil(IndexC./1e7));
IndexP=IndexP+(1-ceil(IndexP./1e7));

IndexC(:,2)=AreasC(:,2);
IndexP(:,2)=AreasP(:,2);
Cimage=zeros(size(Image));
SizeIC=size(IndexC);
SizeIC=SizeIC(1);
SizeIP=size(IndexP);
SizeIP=SizeIP(1);
for i=1:SizeIC
    Cimage(CentralLabel==IndexC(i,2))=IndexC(i,1);
end
for i=1:SizeIP
    Cimage(RingLabel==IndexP(i,2))=IndexP(i,1);
end
max(max(Cimage)); % Flag Alex, scaling of the heatmap
max(IndexC(:,1)); % Flag Alex, scaling of the heatmap
max(IndexP(:,1)); % Flag Alex, scaling of the heatmap
Inversion=zeros(size(Image));
Inversion(Cimage==0)=1;
Counter=1-Inversion;
Cimage=Cimage+Counter.*255;
Image=single(Image).*single(Inversion);
Image=Image+Cimage;
CmapFull=[Cmap1; Cmap2];
axes(handles.axes1); hold off
imshow(Image,CmapFull); hold on
%colorbar % Flag Alex, heatmap colorbar
RGBim=ind2rgb(Image,CmapFull);
set(handles.HoleSizeHeatmap,'UserData',RGBim);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in SquareRootMode.
function SquareRootMode_Callback(hObject, eventdata, handles)
% hObject    handle to SquareRootMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mode=get(handles.SquareRootMode,'Value');
if Mode==0
    set(handles.LinearMode,'Value',1);
    set(handles.LogMode,'Value',0);
    set(handles.SquareRootMode,'UserData',0);
end
if Mode==1
    set(handles.LinearMode,'Value',0);
    set(handles.LogMode,'Value',0);
    set(handles.SquareRootMode,'UserData',1);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in LinearMode.
function LinearMode_Callback(hObject, eventdata, handles)
% hObject    handle to LinearMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mode=get(handles.LinearMode,'Value');
if Mode==0
    set(handles.SquareRootMode,'Value',0);
    set(handles.LogMode,'Value',1);
    set(handles.SquareRootMode,'UserData',2);
end
if Mode==1
    set(handles.SquareRootMode,'Value',0);
    set(handles.LogMode,'Value',0);
    set(handles.SquareRootMode,'UserData',0);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in LogMode.
function LogMode_Callback(hObject, eventdata, handles)
% hObject    handle to LogMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mode=get(handles.LogMode,'Value');
if Mode==0
    set(handles.SquareRootMode,'Value',1);
    set(handles.LogMode,'Value',0);
    set(handles.SquareRootMode,'UserData',1);
end
if Mode==1
    set(handles.SquareRootMode,'Value',0);
    set(handles.LinearMode,'Value',0);
    set(handles.SquareRootMode,'UserData',2);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in SaveImage.
function SaveImage_Callback(hObject, eventdata, handles)
% hObject    handle to SaveImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
OpenData=get(handles.OpenImage,'UserData');
File=OpenData{2};
File=[File '_hole_heatmap.tif'];
Image=get(handles.HoleSizeHeatmap,'UserData');
imwrite(Image,File,'tif');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in ViewPenetration.
function ViewPenetration_Callback(hObject, eventdata, handles)
% hObject    handle to ViewPenetration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Image=get(handles.SelectRectangle,'UserData');

ImageSize=size(Image);
Cmap1=colormap(gray(256));
CutOffDistanceTransforms=get(handles.ViewPenetration,'UserData');
Cimage=CutOffDistanceTransforms{3};
Cmap2=colormap(jet(255));

Inversion=zeros(size(Image));
Inversion(Cimage==0)=1;
Counter=1-Inversion;

Cimage=Cimage+Counter.*255;
Image=single(Image).*single(Inversion);
Image=Image+Cimage;
CmapFull=[Cmap1; Cmap2];
axes(handles.axes1)
imshow(Image,CmapFull);
RGBim=ind2rgb(Image,CmapFull);
CutOffDistanceTransforms{4}=RGBim;
save([tempdir,'CutOff.mat'],'CutOffDistanceTransforms')
set(handles.ViewPenetration,'UserData',CutOffDistanceTransforms);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in BlobCheckbox.
function BlobCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to BlobCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Checking=get(handles.BlobCheckbox,'Value');
if Checking==1
    set(handles.MeshCheckbox,'Value',0);
    set(handles.AutomaticMode,'Value',0);
    set(handles.BlobCheckbox,'UserData',1);
else
    set(handles.MeshCheckbox,'Value',1);
    set(handles.AutomaticMode,'Value',0);
    set(handles.BlobCheckbox,'UserData',0);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in MeshCheckbox.
function MeshCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to MeshCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Checking=get(handles.MeshCheckbox,'Value');
if Checking==1
    set(handles.BlobCheckbox,'Value',0);
    set(handles.AutomaticMode,'Value',0);
    set(handles.BlobCheckbox,'UserData',0);
else
    set(handles.BlobCheckbox,'Value',1);
    set(handles.AutomaticMode,'Value',0);
    set(handles.BlobCheckbox,'UserData',1);
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in SavePenetration.
function SavePenetration_Callback(hObject, eventdata, handles)
% hObject    handle to SavePenetration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
OpenData=get(handles.OpenImage,'UserData');
File=OpenData{2};
File=[File '_penetration_heatmap.tif'];
CutOffDistanceTransforms=get(handles.ViewPenetration,'UserData');
Image=CutOffDistanceTransforms{4};
imwrite(Image,File,'tif');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in ViewBackground.
function ViewBackground_Callback(hObject, eventdata, handles)
% hObject    handle to ViewBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Image=get(handles.SelectRectangle,'UserData');
ImageSize=size(Image);
Mask=get(handles.ViewBackground,'UserData');
axes(handles.axes1); hold off;
NewImage=single(Image).*Mask;
imshow(NewImage,[]);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in RemoveRegion.
function RemoveRegion_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Image=get(handles.SelectRectangle,'UserData');
ImageSize=size(Image);
Mask=get(handles.ViewBackground,'UserData');
axes(handles.axes1); hold off;
NewImage=single(Image).*Mask;
imshow(NewImage,[]);
% Removal=roipoly;
% Mask=Mask-Mask.*Removal;
NewImage=single(Image).*Mask;
imshow(NewImage,[]);
set(handles.ViewBackground,'UserData',Mask)

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in ExecuteAutoThreshold.
function ExecuteAutoThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to ExecuteAutoThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Image=get(handles.SelectRectangle,'UserData');
Background=get(handles.ViewBackground,'UserData');
CellP=handles.ChooseCellPerimeter;
CellMask=CellP{1};
ImgS=smoothn(Image);
axes(handles.axes1); hold off;
imshow(ImgS,[]);
MeanBackgroundIntensity=sum(sum(double(Image).*double(Background)))./sum(sum(Background));
BackgroundValues=double(ImgS(Background==1));
BGmed=median(BackgroundValues);
Subs=min(MeanBackgroundIntensity,BGmed);
ImgS=ImgS-Subs;
ImgS=0.5*(ImgS+abs(ImgS));
Min5=ordfilt2(ImgS,1,ones(5,5));
Filter=fspecial('gaussian',15,2);
MinB=imfilter(Min5,Filter);
NewImage=ImgS-MinB;
NewImage=0.5*(NewImage+abs(NewImage));
NewImage=mat2gray(NewImage);
Xf=fspecial('gaussian',7,1.1);
Apu=deconvlucy(NewImage,Xf);
BackgroundValues=double(Apu(Background==1));
ApuBGmean=sum(sum(Background.*Apu))./sum(sum(Background));
ApuBGmstd=ApuBGmean+2*std(BackgroundValues);
ApuBG67=prctile(BackgroundValues,90);
Subs=max(ApuBGmstd,ApuBG67);
Apu=0.5*((Apu-Subs)+abs(Apu-Subs));
BW=im2bw(Apu,0);
BW=bwlabel(BW);
Stats=regionprops(BW,'Area');
idx=find([Stats.Area]>1);
Automatic=ismember(BW,idx);
Automatic=Automatic.*CellMask;
axes(handles.axes1); hold off;
imshow(Automatic,[])
set(handles.ExecuteAutoThreshold,'UserData',Automatic);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in AutomaticMode.
function AutomaticMode_Callback(hObject, eventdata, handles)
% hObject    handle to AutomaticMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Checking=get(handles.AutomaticMode,'Value');
if Checking==1
    set(handles.BlobCheckbox,'Value',0);
    set(handles.MeshCheckbox,'Value',0);
    set(handles.BlobCheckbox,'UserData',2);
else
    set(handles.BlobCheckbox,'Value',1);
    set(handles.MeshCheckbox,'Value',0);
    set(handles.BlobCheckbox,'UserData',1);
end

Centrals=handles.ChooseCentralZone;
CentralMask=Centrals{1};
RingMask=handles.StoreForRingMask;
size(RingMask);
Positive=get(handles.ExecuteAutoThreshold,'UserData');
Pos1=Positive.*CentralMask;
Pos2=Positive.*RingMask;
set(handles.ThresholdMultiplierValue,'UserData',Pos1);
set(handles.ThresholdSlider,'UserData',Pos2);

% Update handles structure
guidata(hObject, handles);
