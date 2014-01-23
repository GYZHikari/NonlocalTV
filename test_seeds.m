function varargout = test_seeds(varargin)
% TEST_SEEDS MATLAB code for test_seeds.fig
%      TEST_SEEDS, by itself, creates a new TEST_SEEDS or raises the existing
%      singleton*.
%
%      H = TEST_SEEDS returns the handle to a new TEST_SEEDS or the handle to
%      the existing singleton*.
%
%      TEST_SEEDS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEST_SEEDS.M with the given input arguments.
%
%      TEST_SEEDS('Property','Value',...) creates a new TEST_SEEDS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before test_seeds_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test_seeds_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help test_seeds
% GuangyuZhongHikari@gmail.com
% Last Modified by GUIDE v2.5 20-Jan-2014 00:40:02

% Begin initialization code - DO NOT EDIT
addpath(genpath('.'));
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @test_seeds_OpeningFcn, ...
    'gui_OutputFcn',  @test_seeds_OutputFcn, ...
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


% --- Executes just before test_seeds is made visible.
function test_seeds_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to test_seeds (see VARARGIN)

% Choose default command line output for test_seeds
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes test_seeds wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = test_seeds_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in confirm.
function confirm_Callback(hObject, eventdata, handles)
% hObject    handle to confirm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maxSpNum = str2num(get(handles.spnum,'String'));
guidata(hObject,handles);

% --- Executes on button press in openfile
function openfile_Callback(hObject, eventdata, handles)
% hObject    handle to openfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.imdir = '.\img\input\';
handles.spdir = '.\img\output\superpixels\';
handles.saldir = '.\img\output\saliency\';

[handles.imname, handles.imdir] = uigetfile({'*.bmp', 'Image Files'}, 'Pick a file', handles.imdir);
[handles.img] = imread([handles.imdir handles.imname]);
[handles.m,handles.n,handles.k] = size(handles.img);
[handles.spImg, handles.spAdjcMat, handles.spInds, handles.spCnt, handles.spNpx] = gene_superpixel(handles.imdir, handles.imname,handles.maxSpNum, handles.spdir, handles.m, handles.n);
handles.spNum = length(handles.spInds);
handles.sp_im = draw_seed(handles.img, handles.spImg, [],[255,0,0]);
axes(handles.axes1);
imshow(handles.img);
axes(handles.axes1);
imshow(handles.sp_im);
handles.Xrange=get(handles.axes1,'XLim');
handles.Yrange=get(handles.axes1,'YLim');

guidata(hObject,handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.flag = 1;
currPt = get(handles.axes1,'CurrentPoint');
handles.x = currPt(1,1);
handles.y = currPt(1,2);
switch lower(handles.graphMode)
    case 'fore'
        handles.foreseed = [handles.foreseed;currPt(1,1:2)];
    case 'back'
        handles.backseed = [handles.backseed;currPt(1,1:2)];
    case 'eraser'
        handles.color(1,1) = handles.img(round(handles.y),round(handles.x),1);
        handles.color(1,2) = handles.img(round(handles.y),round(handles.x),2);
        handles.color(1,3) = handles.img(round(handles.y),round(handles.x),3);
        handles.color = double(handles.color)./255;
        handles.eraserseed = [handles.eraserseed;currPt(1,1:2)];
        
end
line(handles.x,handles.y,'linewidth', 4,'Color',handles.color);
handles.x0 = handles.x;
handles.y0 = handles.y;
guidata(hObject,handles);


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.flag
    currPt = get(handles.axes1,'CurrentPoint');
    handles.x = currPt(1,1);
    handles.y = currPt(1,2);
    %     handles.Position = [handles.Position; currPt];
    switch lower(handles.graphMode)
        case 'fore'
            handles.foreseed = [handles.foreseed;currPt(1,1:2)];
        case 'back'
            handles.backseed = [handles.backseed;currPt(1,1:2)];
        case 'eraser'
            handles.color(1,1) = handles.img(round(handles.y),round(handles.x),1);
            handles.color(1,2) = handles.img(round(handles.y),round(handles.x),2);
            handles.color(1,3) = handles.img(round(handles.y),round(handles.x),3);
            handles.color = double(handles.color)./255;
            handles.eraserseed = [handles.eraserseed;currPt(1,1:2)];
            
    end
    line([handles.x0 handles.x],[handles.y0 handles.y],'linewidth', 4,'Color',handles.color);
    handles.x0 = handles.x;
    handles.y0 = handles.y;
end
guidata(hObject,handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.flag = 0;
if ~isempty(handles.spImg)
    
    switch lower(handles.graphMode)
        case 'fore'
            handles.foreseed(handles.foreseed<1) = 1; % y in the img
            
            if ~isempty(handles.foreseed)
                handles.foreseed(handles.foreseed(:,2)>handles.m) = handles.m; % y in the img
                handles.foreseed(handles.foreseed(:,1)>handles.n) = handles.n; % x in the img
            end
            handles.foreseed = fix(handles.foreseed);
            index=sub2ind(size(handles.spImg),handles.foreseed(:,2),handles.foreseed(:,1));
            handles.ForeSeed = unique(handles.spImg(index));
            
            commen = intersect(handles.BackSeed,handles.ForeSeed);
            handles.BackSeed = setdiff(handles.BackSeed,commen);
            ind = [];
            if ~isempty(commen)
                backindex = sub2ind(size(handles.spImg),handles.backseed(:,2),handles.backseed(:,1));
                for i = 1:length(commen)
                    [in] = find(handles.spImg(backindex)==commen(i));
                    [ind] = [ind;in];
                end
                handles.backseed(ind,:) = [];
                
            end
            
            sp_im = draw_seed(handles.img, handles.spImg, handles.ForeSeed,[255,0,0]);
            sp_im = draw_seed(sp_im, handles.spImg, handles.BackSeed,[0,0,255]);
            handles.sp_im = sp_im;
            axes(handles.axes1);
            imshow(handles.sp_im);
        case 'back'
            handles.backseed(handles.backseed<1) = 1; % y in the img
            
            if ~isempty(handles.backseed)
                handles.backseed(handles.backseed(:,2)>handles.m) = handles.m; % y in the img
                handles.backseed(handles.backseed(:,1)>handles.n) = handles.n; % x in the img
            end
            handles.backseed = fix(handles.backseed);
            
            index=sub2ind(size(handles.spImg),handles.backseed(:,2),handles.backseed(:,1));
            handles.BackSeed = unique(handles.spImg(index));
            
            commen = intersect(handles.BackSeed,handles.ForeSeed);
            handles.ForeSeed = setdiff(handles.ForeSeed,commen);
            ind = [];
            if ~isempty(commen)
                foreindex = sub2ind(size(handles.spImg),handles.foreseed(:,2),handles.foreseed(:,1));
                for i = 1:length(commen)
                    [in] = find(handles.spImg(foreindex)==commen(i));
                    ind = [ind;in];
                end
                handles.foreseed(ind,:) = [];
            end
            sp_im = draw_seed(handles.img, handles.spImg, handles.ForeSeed,[255,0,0]);
            sp_im = draw_seed(sp_im, handles.spImg, handles.BackSeed,[0,0,255]);
            handles.sp_im = sp_im;
            axes(handles.axes1);
            imshow(handles.sp_im);
            
        case 'eraser'
            handles.eraserseed(handles.eraserseed<1) = 1; % y in the img
            
            if ~isempty(handles.eraserseed)
                handles.eraserseed(handles.eraserseed(:,2)>handles.m) = handles.m; % y in the img
                handles.eraserseed(handles.eraserseed(:,1)>handles.n) = handles.n; % x in the img
            end
            handles.eraserseed = fix(handles.eraserseed);
            
            index=sub2ind(size(handles.spImg),handles.eraserseed(:,2),handles.eraserseed(:,1));
            handles.EraSeed = unique(handles.spImg(index));
            
            commen = intersect(handles.ForeSeed,handles.EraSeed);
            handles.ForeSeed = setdiff(handles.ForeSeed,commen);
            ind = [];
            if ~isempty(commen)
                foreindex = sub2ind(size(handles.spImg),handles.foreseed(:,2),handles.foreseed(:,1));
                for i = 1:length(commen)
                    [in] = find(handles.spImg(foreindex)==commen(i));
                    ind = [ind;in];
                end
                handles.foreseed(ind,:) = [];
            end
            
            commen = intersect(handles.BackSeed,handles.EraSeed);
            handles.BackSeed = setdiff(handles.BackSeed,commen);
            ind = [];
            if ~isempty(commen)
                backindex = sub2ind(size(handles.spImg),handles.backseed(:,2),handles.backseed(:,1));
                for i = 1:length(commen)
                    [in] = find(handles.spImg(backindex)==commen(i));
                    [ind] = [ind;in];
                end
                handles.backseed(ind,:) = [];
                
            end
            
            
            sp_im = draw_seed(handles.img, handles.spImg, handles.ForeSeed,[255,0,0]);
            sp_im = draw_seed(sp_im, handles.spImg, handles.BackSeed,[0,0,255]);
            handles.sp_im = sp_im;
            axes(handles.axes1);
            imshow(handles.sp_im);
            handles.EraSeed = [];
            handles.eraserseed =[];
            
    end
end
guidata(hObject,handles);


% --- Executes on selection change in brush.
function brush_Callback(hObject, eventdata, handles)
% hObject    handle to brush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns brush contents as cell array
%        contents{get(hObject,'Value')} returns selected item from brush
str = get(handles.brush,'string');
index = get(handles.brush,'value');
str1 = char(str(index));
switch(str1)
    case 'ForeSeed'
        handles.graphMode = 'fore';
        handles.color = [1,0,0];
    case 'BackSeed'
        handles.graphMode = 'back';
        handles.color = [0,0,1];
    case 'Eraser'
        handles.graphMode = 'eraser';
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function brush_CreateFcn(hObject, eventdata, handles)
% hObject    handle to brush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
handles.flag = 0;
handles.maxSpNum = 200;
handles.spImg = [];
handles.foreseed= [];
handles.backseed =[];
handles.BackSeed = [];
handles.ForeSeed = [];
handles.EraSeed = [];
handles.eraserseed =[];
handles.graphMode = 'fore';
handles.color = [1,0,0];
guidata(hObject,handles);

% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in solver.
function solver_Callback(hObject, eventdata, handles)
% hObject    handle to solver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = get(handles.solver,'string');
index = get(handles.solver,'value');
str1 = char(str(index));
switch(str1)
    case 'Manifold ranking'
        handles.solver = 'manirank';
    case 'nonlocal_TV_h1_iter'
        handles.solver = 'H1';
    case 'nonlocal_TV_ADM'
        handles.solver = 'ADM';
end
guidata(hObject,handles);

% Hints: contents = cellstr(get(hObject,'String')) returns solver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from solver


% --- Executes on button press in compute.
function compute_Callback(hObject, eventdata, handles)
% hObject    handle to compute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global param;
param.featMode = 1;
param.connectMode = 'local';
param.guass = 0;
param.featDistOpt = 'MVSSER';
param.simiMode = 'similar';
param.simiParam = -10;

handles.prior = zeros(handles.spNum,1);
handles.prior(handles.ForeSeed) = 1;

[handles.labFea,~] = gene_feature(handles.img, handles.spImg, handles.spCnt, handles.spNpx, param);
handles.labFea = normalize(handles.labFea);
handles.edges = gene_connect_edges(handles.spAdjcMat,handles.spNum,handles.BackSeed,[],2,param,handles.labFea,handles.spCnt);
handles.affmat = gene_affmat(handles.edges,handles.spNum,handles.spCnt,handles.labFea,param);
handles.W = reshape(handles.affmat',1,handles.spNum.^2);

handles.sal = solve_NLH1_iter(handles.W,handles.prior',handles.spNum,10,2,1,1,1);
handles.salMap = saliencySp2img(handles.sal,handles.spImg);
axes(handles.axes2);
imshow(handles.salMap);
imwrite(handles.salMap,[handles.saldir,handles.imname]);