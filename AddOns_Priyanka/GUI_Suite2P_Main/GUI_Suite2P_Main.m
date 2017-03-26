function varargout = GUI_Suite2P_Main(varargin)
% GUI_SUITE2P_MAIN MATLAB code for GUI_Suite2P_Main.fig
%      GUI_SUITE2P_MAIN, by itself, creates a new GUI_SUITE2P_MAIN or raises the existing
%      singleton*.
%
%      H = GUI_SUITE2P_MAIN returns the handle to a new GUI_SUITE2P_MAIN or the handle to
%      the existing singleton*.
%
%      GUI_SUITE2P_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SUITE2P_MAIN.M with the given input arguments.
%
%      GUI_SUITE2P_MAIN('Property','Value',...) creates a new GUI_SUITE2P_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Suite2P_Main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Suite2P_Main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Suite2P_Main

% Last Modified by GUIDE v2.5 02-Feb-2017 15:53:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Suite2P_Main_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Suite2P_Main_OutputFcn, ...
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


% --- Executes just before GUI_Suite2P_Main is made visible.
function GUI_Suite2P_Main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Suite2P_Main (see VARARGIN)

% Choose default command line output for GUI_Suite2P_Main
handles.output = hObject;

% get computer name
!hostname > hostname.txt
handles.computername = textread('hostname.txt','%s');

% load computer specific settings
handles = My_FilePaths(handles);
if ~isempty(handles.ops0.toolbox_path)
    % toolbox
    if exist(handles.ops0.toolbox_path, 'dir')
        addpath(handles.ops0.toolbox_path) % add local path to the toolbox
    else
        error('toolbox_path does not exist, please change toolbox_path');
    end
end
% reposition GUI
movegui(hObject,'northwest'); 

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Suite2P_Main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Suite2P_Main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SetUpSession_CLEAR.
function SetUpSession_CLEAR_Callback(hObject, eventdata, handles)
% hObject    handle to SetUpSession_CLEAR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SetUpSession_TotalSessions.Data = [0 1 1]';
handles.db = [];
handles.session_list.String = {' '};
guidata(hObject, handles);


% --- Executes on button press in SetUpSession_AddNew.
function SetUpSession_AddNew_Callback(hObject, eventdata, handles)
% hObject    handle to SetUpSession_AddNew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = uigetfile_n_dir(char(handles.filepaths.Data(1)),'Select experimental session folder');

if ~isempty (X)
    handles.SetUpSession_TotalSessions.Data(1) = handles.SetUpSession_TotalSessions.Data(1) + 1;
    y = handles.SetUpSession_TotalSessions.Data(1);
    M =  regexp(char(X(1)),filesep,'split');
    handles.db(y).mouse_name    = char(M(end-2));
    handles.db(y).date = char(M(end-1));
    %handles.db(y).expts =  str2num(char(M(end)));
    handles.db(y).expts(1) =  M(end);
    handles.db(y).nplanes = handles.SetUpSession_TotalSessions.Data(2);
    % for red channel
%     handles.db(y).AlignToRedChannel = handles.AlignToRed.Value;
%     handles.db(y).nchannels_red = 2; % how many channels did the red block have in total (assumes red is last)
    handles.session_list.String(y:y+size(X,2)-1) = X;
    if size(X,2)>1
        for i = 2:size(X,2)
            M =  regexp(char(X(i)),filesep,'split');
            %handles.db(y).expts = [handles.db(y).expts str2num(char(M(end)))];
            handles.db(y).expts(i) = M(end);
%             if handles.AlignToRed.Value
%                 handles.db(y).expred(i) = M(end);
%             end
        end
    end
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in make_master_file.
function make_master_file_Callback(hObject, eventdata, handles)

%% root paths for files and temporary storage (ideally an SSD drive)
% Suite2P assumes a folder structure, check out README file
handles.ops0.RootStorage = char(handles.filepaths.Data(1));
% copy each remote tiff locally first, into this file
handles.ops0.temp_tiff = fullfile(char(handles.filepaths.Data(2)),'temp.tif');
% location for binary file
handles.ops0.RegFileRoot = char(handles.filepaths.Data(2)); 
% set to 1 for batch processing on a limited hard drive
handles.ops0.DeleteBin = handles.limited_hard_drive.Value;
% a folder structure is created inside
handles.ops0.ResultsSavePath = char(handles.filepaths.Data(3)); 
% leave empty to NOT save registered tiffs (slow)
if handles.save_registered_tiffs.Value
    handles.ops0.RegFileTiffLocation = char(handles.filepaths.Data(4));
else
    handles.ops0.RegFileTiffLocation = char({''});
end

%% registration settings
% if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. 
% You only need the Nvidia drivers installed (not CUDA).
handles.ops0.useGPU = handles.use_gpu.Value;
% skip (0) if data is already registered
handles.ops0.doRegistration = handles.register_data.Value;
% shows the image targets for all planes to be registered
handles.ops0.showTargetRegistration = handles.show_reference_image.Value;
% set to 0 for non-whitened cross-correlation
handles.ops0.PhaseCorrelation = handles.non_whitened_cross_correlation.Value;
% sub-pixel alignment
if  handles.sub_pixel_align.Value
    handles.ops0.SubPixel = Inf; % Inf is the exact number from phase correlation
else
    handles.ops0.SubPixel = 2; % 2 is alignment by 0.5 pixel
end
% number of images to include in the first registration pass
handles.ops0.NimgFirstRegistration  =  handles.registration_settings.Data(1); 
% upsampling factor during registration, 
% 1 for no upsampling is much faster, 2 may give better subpixel accuracy
handles.ops0.registrationUpsample   = handles.upsample.Value + 1;

% use red channel
handles.ops0.AlignToRedChannel = handles.AlignToRed.Value;
handles.ops0.nchannels = handles.SetUpSession_TotalSessions.Data(3);

% % currently unused
% handles.ops0.nimgbegend = handles.registration_settings_unused.Data(1); % how many frames to average at the beginning and end of each experiment
% handles.ops0.NiterPrealign = handles.registration_settings_unused.Data(2);

%% Roi detection settings
handles.ops0.getROIs = handles.getROIs.Value;
handles.ops0.clustModel = handles.Cluster_Mode.SelectedObject.String; % standard or neuropil
handles.ops0.neuropilSub = handles.Neuropil_Subtraction_Model.SelectedObject.String; % none, surround or model
% during optimization, show a figure of the clusters
handles.ops0.ShowCellMap = handles.show_cell_map.Value;
% how many clusters to start with
handles.ops0.Nk0 = handles.cluster_settings.Data(1);
% how many clusters to end with (before anatomical segmentation)
handles.ops0.Nk = handles.cluster_settings.Data(2);
% spatial smoothing length for clustering; encourages localized clusters
handles.ops0.sig = handles.cluster_settings.Data(3);
% how many SVD components for cell clustering
handles.ops0.nSVDforROI = handles.cluster_settings.Data(4);
% how many (binned) timepoints to do the SVD based on
handles.ops0.NavgFramesSVD = handles.cluster_settings.Data(5); 
% expected diameter of cells (used for 0.25 * pi/4*diam^2 < npixels < 10*pi/4*diam^2)
handles.clustrules.diameter = handles.cluster_settings.Data(6);

% % currently unused
% handles.ops0.niterclustering = handles.cluster_settings.Data(4); % how many iterations of clustering
% handles.ops0.getSVDcomps = handles.get_svd.Value;
% handles.ops0.nSVD = handles.cluster_settings.Data(7); % how many SVD components to keep
% % these are modifiable settings for classifying ROIs post-clustering, these settings can be over-ridden in the GUI after running the pipeline
% handles.ops0.LoadRegMean = handles.load_registration_mean.Value;
% handles.clustrules.MaxNpix = handles.deconvolution_settings.Data(1);
% handles.clustrules.MinNpix = handles.deconvolution_settings.Data(2);
% handles.clustrules.Compact = handles.deconvolution_settings.Data(3);
% handles.clustrules.parent.minPixRelVar = handles.deconvolution_settings.Data(4);
% handles.clustrules.parent.PixelFractionThreshold = handles.deconvolution_settings.Data(5);
% handles.clustrules.parent.MaxRegions = handles.deconvolution_settings.Data(6);

%% Spike deconvolution options
% imaging rate (cumulative over planes!). 
% Approximate, for initialization of deconvolution kernel.
handles.ops0.imageRate = handles.deconvolution_settings.Data(1);
% decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.
handles.ops0.sensorTau = handles.deconvolution_settings.Data(2);
% for the neuropil contamination to be less than this (sometimes good, i.e. for interneurons)
handles.ops0.maxNeurop = handles.deconvolution_settings.Data(3);
% whether to re-estimate kernel during optimization
% (default kernel is "reasonable", if you give good timescales)
handles.ops0.recomputeKernel = handles.recompute_kernel.Value;
% whether the same kernel should be estimated for all neurons
% (robust, only set to 0 if SNR is high and recordings are long)
handles.ops0.sameKernel = handles.same_kernel.Value;

handles.db0 = handles.db;
%handles.make_master_file.String = 'running ...';
%handles.make_master_file.BackgroundColor = [0.5 0.9400 0.9400];
set(hObject,'BackgroundColor','cyan','String','running ...');
pause(0.5);
guidata(hObject, handles);

%% RUN THE PIPELINE HERE
for iexp = 1:length(handles.db)   
    handles.session_list.Value = iexp;
    display('running pipeline....');
    run_pipeline(handles.db(iexp), handles.ops0, handles.clustrules);
    
    if handles.deconvolve_spikes.Value
        display('deconvolving spikes....');
        set(hObject,'String','deconvolving ...');
        pause(0.5);
        % deconvolved data into (dat.)cl.dcell, and neuropil subtraction coef
        % MAKE SURE YOU COMPILE THIS FIRST FOR DECONVOLUTION
        % mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp)
        add_deconvolution(handles.ops0, handles.db0(iexp), handles.clustrules);
    end
    
    % add red channel information (if it exists)
     % run_REDaddon(iexp, db, ops0);
end
display('done!');
% handles.make_master_file.String = 'GO';
% handles.make_master_file.BackgroundColor = [0.9400 0.9400 0.9400];
set(hObject,'BackgroundColor',[0.9400 0.9400 0.9400],'String','GO');
pause(0.5);
guidata(hObject, handles);


% --- Executes on button press in AlignToRed.
function AlignToRed_Callback(hObject, eventdata, handles)
% hObject    handle to AlignToRed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AlignToRed
