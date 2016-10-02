function varargout = ROI_prune(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ROI_prune_OpeningFcn, ...
    'gui_OutputFcn',  @ROI_prune_OutputFcn, ...
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

% --- Executes just before ROI_prune is made visible.
function ROI_prune_OpeningFcn(hObject, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% varargin   command line arguments to ROI_prune (see VARARGIN)

h.output = hObject;
h.ROIs_compare = [];
h.polygon = [];
my_axes(1) = h.axes4; % F trace
my_axes(2) = h.axes6; % spike times
linkaxes(my_axes,'x');
guidata(hObject, h);
whitebg(gcf,'k');
set(gcf,'color','black');

% reposition GUI
movegui(hObject,'northwest'); 
% UIWAIT makes ROI_prune wait for user response (see UIRESUME)
% uiwait(h.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ROI_prune_OutputFcn(hObject, eventdata, h)
% varargout  cell array for returning output args (see VARARGOUT);
varargout{1} = h.output;

function Load_plane_Callback(hObject, eventdata, h)
flag = 0;
try
    if isfield(h, 'dat') && isfield(h.dat, 'filename')
        root = fileparts(h.dat.filename);
    else
        root = 'D:\DATA\F\';
    end
    [filename1,filepath1]=uigetfile(root, 'Select Data File');
    h.dat = load(fullfile(filepath1, filename1));
    % if its a processed session file - extract session params
    % and load the original file
    %     if isfield(h.dat,'S')
    %         sessionparams = h.dat.S;
    %         h = rmfield(h,'dat');
    %         h.dat.dat = load(sessionparams.filename);
    %         h.dat.dat.cl = sessionparams.cl;
    %         h.dat.dat.res = sessionparams.res;
    %     end
    
    set(h.figure1, 'Name', filename1);
    flag = 1;
catch
end

if flag
    % if the user selected a file, do all the initializations
    rng('default')
    
    newsession = 1;
    if isfield(h.dat, 'dat')
        session = h.dat.dat;
        h = rmfield(h,'dat');
        
        % first load the master file
        h.dat = load(session.filename);
        
        % copy over fields
        sessionfields = fieldnames(session);
        for i = 1:size(sessionfields,1)
            if isfield(h.dat,sessionfields(i))
                h.dat.(sessionfields{i}) = session.(sessionfields{i});
            end
        end
        h.dat.F.ichosen = 1;
        
        clear session sessionfields
        
        %h.dat = h.dat.dat;
        %         h = splitROIleftright_PG(h);
        %         h = buildHue(h);
        %         h = buildLambdaValue(h);
        %         h = update_cell_neuropil_traces(h);
        newsession = 0;
    end
    
    
    if newsession
        h.dat.filename = fullfile(filepath1, filename1);
        h.dat.cl.Mrs      = [h.dat.stat.mrs]./[h.dat.stat.mrs0];
        h.dat.cl.npix     = [h.dat.stat.npix];
        h.dat.cl.Ly       = numel(h.dat.ops.yrange);
        h.dat.cl.Lx       = numel(h.dat.ops.xrange);
        h.dat.cl.MeanM    = 2*mean(h.dat.res.M);
        h.dat.cl.excluded_pixels  = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
        h.dat.cl.excluded_regions = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
        h.dat.cl.excl_pix_perc    = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
        h.dat.cl.topregion        = ones(h.dat.cl.Ly, h.dat.cl.Lx);
        h.dat.res.iclust = reshape(h.dat.res.iclust, h.dat.cl.Ly, h.dat.cl.Lx);
        Nk = h.dat.ops.Nk;
        h.dat.ops.Nk = numel(h.dat.stat);
        h.dat.cl.rands_orig   = .1 + .8 * rand(1, h.dat.ops.Nk);
        h.dat.cl.rands        = h.dat.cl.rands_orig;
        h = get_parent_stats(h);

        if isfield(h.dat, 'clustrules')
            % ROI rules
            h.dat.res.Mrs_thresh_orig = h.dat.clustrules.Compact;
            h.dat.cl.npix_low_orig    = h.dat.clustrules.MinNpix;
            h.dat.cl.npix_high_orig   = h.dat.clustrules.MaxNpix;
        else
            % ROI rules
            h.dat.res.Mrs_thresh_orig   = 3;
            h.dat.cl.npix_low_orig      = 20;
            h.dat.cl.npix_high_orig     = 500;
        end
        
        % parent rules
        h.dat.cl.mrs_parent_max = Inf;
        h.dat.cl.npix_res_max   = Inf;
        h.dat.cl.npix_par_max   = Inf;
        h.dat.cl.nreg_max       = Inf;
        h.dat.cl.VperPix_min    = 0;
        
        h = setOriginalThresh(h);
    end
    
    % start with unit vector map
    lam = h.dat.res.lambda;
    h.dat.img0.V = max(0, min(1, .5 * reshape(lam, h.dat.cl.Ly, h.dat.cl.Lx)/mean(lam(:))));
    
    h.dat.ylim = [0 h.dat.cl.Ly];
    h.dat.xlim = [0 h.dat.cl.Lx];
    
    if newsession
        h.dat.cl.manual  = zeros(h.dat.ops.Nk, 1);
        h.dat.cl.redcell = zeros(h.dat.ops.Nk, 1);
    end
    h                = splitROIleftright_PG(h);
    
    icell = find(h.dat.cl.iscell);
    if ~isempty(icell)
        h.dat.F.ichosen = icell(1); %ceil(rand * numel(icell))
    else
        h.dat.F.ichosen = 1; %ceil(rand * numel(icell))
    end
    
    Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
    Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
    h.dat.img1.Sat     = Sat;
    h.dat.img2.Sat     = Sat;
    
    h = buildHue(h);
    h = buildLambdaValue(h);
    
    % x and y limits on subquadrants
    h.dat.figure.x0all = round(linspace(0, 19/20*h.dat.cl.Lx, 4));
    h.dat.figure.y0all = round(linspace(0, 19/20*h.dat.cl.Ly, 4));
    h.dat.figure.x1all = round(linspace(1/20 * h.dat.cl.Lx, h.dat.cl.Lx, 4));
    h.dat.figure.y1all = round(linspace(1/20 * h.dat.cl.Ly, h.dat.cl.Ly, 4));
    
    h = update_cell_neuropil_traces(h);
end

% update parameters display on the GUI
set(h.parent_selection_compactness,'String',num2str(h.dat.cl.mrs_parent_max));
set(h.parent_selection_max_pixel_count,'String',num2str(h.dat.cl.npix_par_max));
set(h.parent_selection_max_pixel_residual,'String',num2str(h.dat.cl.npix_res_max));
set(h.parent_selection_max_region_count,'String',num2str(h.dat.cl.nreg_max));
set(h.parent_selection_min_pixel_variance,'String',num2str(h.dat.cl.VperPix_min));

% roi inclusion rules - update displayed values
set(h.roi_selection_compactness,'String', num2str(h.dat.res.Mrs_thresh));
set(h.roi_selection_max_pixel_count,'String', num2str(h.dat.cl.npix_high));
set(h.roi_selection_min_pixel_count,'String', num2str(h.dat.cl.npix_low));

% set all quadrants as not visited
h.quadvalue = zeros(3);
for j = 1:3
    for i = 1:3
        set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor',[0 0 0]);
    end
end

h.dat.maxmap = 1;
ops = h.dat.ops;
if isfield(ops, 'mimg1') && ~isempty(ops.mimg1)
    h.dat.maxmap = h.dat.maxmap + 1;
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimg1(ops.yrange, ops.xrange);
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end
if isfield(ops, 'mimgRED') && ~isempty(ops.mimgRED)
    h.dat.maxmap = h.dat.maxmap + 1;
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimgRED(ops.yrange, ops.xrange);
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end

h.dat.procmap = 0;
h.dat.map = 1;

h.window_size.String = num2str(size(h.dat.F.trace,2));
redraw_fluorescence_PG(h);
update_display_mode(hObject, h); %redraw_figure(h);

disp('session loaded');

guidata(hObject,h)

function [h] = update_cell_neuropil_traces(h)
h.dat.F.Fcell = h.dat.Fcell;
% h.dat.Fcell = [];    % PG
if isfield(h.dat, 'FcellNeu')
    h.dat.F.FcellNeu = h.dat.FcellNeu;
    % h.dat.FcellNeu = [];
    %     if mean(sign(h.dat.F.FcellNeu{1}))<0
    %         for j = 1:length(h.dat.F.FcellNeu)
    %             h.dat.F.FcellNeu{j} = - h.dat.F.FcellNeu{j};
    %             h.dat.F.Fcell{j} = h.dat.F.Fcell{j} + h.dat.F.FcellNeu{j};
    %         end
    %     end
    if isfield(h.dat.cl, 'dcell')
        for k = 1:length(h.dat.F.FcellNeu)
            cellcount = 0;
            for j = 1:size(h.dat.FcellNeu{k},1) % no. of ROIs
                if h.dat.cl.isroi(k,j)
                    cellcount = cellcount + 1;
                    if isfield(h.dat.cl.dcell{cellcount}, 'B')
                        if ~h.neuropil_subtracted_c.Value
                            c2 = h.dat.cl.dcell{cellcount}.B(3);
                        else
                            c2 = min(1,h.dat.cl.dcell{cellcount}.B(3));
                        end
                        c1 = h.dat.cl.dcell{cellcount}.B(2);
                        h.dat.F.FcellNeu{k}(j,:) = c1 + c2 * h.dat.F.FcellNeu{k}(j,:);
                    end
                end
            end
        end
    end
    h.dat.F.trace = [];
end
for i = 1:length(h.dat.F.Fcell)
    h.dat.F.trace = cat(2, h.dat.F.trace, h.dat.F.Fcell{i});
end
if isfield(h.dat.F, 'FcellNeu')
    h.dat.F.neurop = [];
    for i = 1:length(h.dat.F.FcellNeu)
        h.dat.F.neurop = cat(2, h.dat.F.neurop, h.dat.F.FcellNeu{i});
    end
end

function Update_MaskType_Callback(hObject, eventdata, h)
switch  h.mask_type.SelectedObject.String
    case 'Unit Vectors'
        % binary mask
        h.dat.img0.V = ones(h.dat.cl.Ly, h.dat.cl.Lx);
    case 'VarianceExp'
        % variance explained mask
        h.dat.img0.V = reshape(h.dat.res.M, h.dat.cl.Ly, h.dat.cl.Lx) / h.dat.cl.MeanM;
    case 'Binary'
        % unit vector mask
        h.dat.img0.V = 10 * reshape(h.dat.res.lambda, h.dat.cl.Ly, h.dat.cl.Lx);
end
h = buildLambdaValue(h);
guidata(hObject,h);
update_display_mode(hObject, h); %redraw_figure(h);

function Update_HueType_Callback(hObject, eventdata, h)
switch  h.mask_type.SelectedObject.String
    case 'Randomize'
        % randomize hue
        rng('shuffle')
        h.dat.cl.rands     = rand(1, h.dat.ops.Nk);
        h.dat.cl.rands(1)  = .15;
        h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
        h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
    case 'Default'
        % original default hue
        h.dat.cl.rands   = h.dat.cl.rands_orig;
        h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
        h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
end
guidata(hObject,h);
update_display_mode(hObject, h); %redraw_figure(h);

function update_ROI_exclusion_Callback(hObject, eventdata, h)
h.dat.res.Mrs_thresh = str2double(get(h.roi_selection_compactness,'String'));
h.dat.cl.npix_high = str2double(get(h.roi_selection_max_pixel_count,'String'));
h.dat.cl.npix_low = str2double(get(h.roi_selection_min_pixel_count,'String'));

h.dat.cl.mrs_parent_max = str2double(get(h.parent_selection_compactness,'String'));
h.dat.cl.npix_par_max = str2double(get(h.parent_selection_max_pixel_count,'String'));
h.dat.cl.npix_res_max = str2double(get(h.parent_selection_max_pixel_residual,'String'));
h.dat.cl.nreg_max = str2double(get(h.parent_selection_max_region_count,'String'));
h.dat.cl.VperPix_min = str2double(get(h.parent_selection_min_pixel_variance,'String'));

h = splitROIleftright_PG(h);

h = buildLambdaValue(h);
update_display_mode(hObject, h); %redraw_figure(h);
guidata(hObject,h);

function h = Manual_selection_Callback(hObject, eventdata, h)
[x, y] = ginputc(1, 'color','w');
x = min(max(1, round(x)), h.dat.cl.Lx);
y = min(max(1, round(y)), h.dat.cl.Ly);
h.dat.F.ichosen = h.dat.res.iclust(y, x);
redraw_fluorescence_PG(h);

Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
h.dat.img1.Sat     = Sat;
h.dat.img2.Sat     = Sat;
h = buildLambdaValue(h);

guidata(hObject,h);
update_display_mode(hObject, h); %redraw_figure(h);

function Save_proc_file_Callback(hObject, eventdata, h)
% h.dat.F.trace = [];
% dat = h.dat;
%[file, path] = uiputfile([h.dat.filename(1:end-4) '_proc.mat'],'Save Session As');
%save(fullfile(path, file), 'dat');
Save_my_traces(h);

function figure1_ResizeFcn(hObject, eventdata, h)

function Quadrant_Callback(hObject, eventdata, h)
iy = str2num(eventdata.Source.String(2));
ix = str2num(eventdata.Source.String(3));
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

function paint_quadbutton(h, iy, ix)
for j = 1:3
    for i = 1:3
        if h.quadvalue(j,i)==1
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor','yellow');
        end
    end
end
set(h.(sprintf('Q%d%d', iy,ix)), 'BackgroundColor','red');

% --- Executes on button press in Full_View.
function Full_View_Callback(hObject, eventdata, h)
h.dat.ylim = [0 h.dat.cl.Ly];
h.dat.xlim = [0 h.dat.cl.Lx];
guidata(hObject,h);
update_display_mode(hObject, h); %redraw_figure(h);

function quadrant(hObject, h, iy, ix)
h.dat.ylim = [h.dat.figure.y0all(iy) h.dat.figure.y1all(iy+1)];
h.dat.xlim = [h.dat.figure.x0all(ix) h.dat.figure.x1all(ix+1)];
h.quadvalue(iy, ix) = 1;
guidata(hObject,h);
update_display_mode(hObject, h); %redraw_figure(h);

% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, h)
switch eventdata.Key
    case 'f'
        % flip currently selected unit
        h.dat.cl.manual(h.dat.F.ichosen) = .5 - h.dat.cl.iscell(h.dat.F.ichosen);
        h = splitROIleftright_PG(h);
        h = buildLambdaValue(h);
        guidata(hObject,h);
        if h.dat.maxmap==1
            update_display_mode(hObject, h); %redraw_figure(h);
        end
    case 's'
        % manual selection of units
        Manual_selection_Callback(hObject, eventdata, h);
    case 'q'
        h.display_mode.Value = 4;
        update_display_mode(hObject, h);
    case 'w'
        h.display_mode.Value = 1;
        update_display_mode(hObject, h);
    case 'e'
        h.dat.map = 3;
        if h.dat.maxmap>2
            h.display_mode.Value = 2;
            update_display_mode(hObject, h);
        end
    case 'p'
        h.display_mode.Value = 3;
        update_display_mode(hObject, h);
    case 'rightarrow'
        roi_list = circshift(h.dat.cl.iscell,size(h.dat.cl.iscell,2)-h.dat.F.ichosen,2);
        roi_id = find(roi_list,1) + h.dat.F.ichosen;
        h.dat.F.ichosen = mod(roi_id,size(h.dat.cl.iscell,2));
        
        Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
        Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
        h.dat.img1.Sat     = Sat;
        h.dat.img2.Sat     = Sat;
        h = buildLambdaValue(h);
        guidata(hObject,h);
        update_display_mode(hObject, h);
        redraw_fluorescence_PG(h);
    case 'leftarrow'
        roi_list = circshift(h.dat.cl.iscell,size(h.dat.cl.iscell,2)-h.dat.F.ichosen,2);
        roi_id =h.dat.F.ichosen - (size(h.dat.cl.iscell,2) -  find(roi_list,2,'last') );
        h.dat.F.ichosen = mod(roi_id(1),size(h.dat.cl.iscell,2));
        
        Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
        Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
        h.dat.img1.Sat     = Sat;
        h.dat.img2.Sat     = Sat;
        h = buildLambdaValue(h);
        guidata(hObject,h);
        update_display_mode(hObject, h);
        redraw_fluorescence_PG(h);
end

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, h)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

z = round(eventdata.Source.CurrentAxes.CurrentPoint(1,:));
x = round(z(1));
y  = round(z(2));

if ~h.define_polygon.Value
    if (x<=h.dat.cl.Lx) && (y<=h.dat.cl.Ly)
        % disp(eventdata.Source.SelectionType)
        % keyboard;
        % manual selection of units
        x = min(max(1, round(x)), h.dat.cl.Lx);
        y = min(max(1, round(y)), h.dat.cl.Ly);
        h.dat.F.ichosen = h.dat.res.iclust(y, x);
        
        switch eventdata.Source.SelectionType
            case 'alt'
                % flip currently selected unit
                h.dat.cl.manual(h.dat.F.ichosen) = .5 - h.dat.cl.iscell(h.dat.F.ichosen);
                h = splitROIleftright_PG(h);
                h = buildLambdaValue(h);
            case 'open'
                % unpin the manual selection on this cell
                h.dat.cl.manual(h.dat.F.ichosen) = 0;
                h = splitROIleftright_PG(h);
                h = buildLambdaValue(h);
            case 'extend'
                h.dat.cl.redcell(h.dat.F.ichosen) = 1 -  h.dat.cl.redcell(h.dat.F.ichosen);
                
                if h.dat.cl.redcell(h.dat.F.ichosen)==1
                    h.dat.cl.rands(h.dat.F.ichosen) = 0;
                else
                    h.dat.cl.rands(h.dat.F.ichosen) = h.dat.cl.rands_orig(h.dat.F.ichosen);
                end
                h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
                h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
                
                if h.dat.cl.redcell(h.dat.F.ichosen)
                    display('red')
                else
                    display('not red')
                end
        end
        
        redraw_fluorescence_PG(h);
        
        % update ROI stats
        %     if h.dat.cl.iscell(h.dat.F.ichosen)
        set(h.current_roi_compactness,'String', num2str(h.dat.cl.Mrs(h.dat.F.ichosen),3));
        set(h.current_roi_pixel_count,'String', num2str(h.dat.cl.npix(h.dat.F.ichosen)));
        set(h.current_roi_pixel_residual,'String', num2str(h.dat.cl.npix_res(h.dat.F.ichosen)));
        set(h.current_roi_region_count,'String', num2str(h.dat.cl.nreg(h.dat.F.ichosen)));
        set(h.current_roi_pixel_variance,'String', num2str(h.dat.cl.VperPix(h.dat.F.ichosen),2));
        %     end
        
        Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
        Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
        h.dat.img1.Sat     = Sat;
        h.dat.img2.Sat     = Sat;
        h = buildLambdaValue(h);
        
        guidata(hObject,h);
        update_display_mode(hObject, h); %redraw_figure(h);
    end
else
    h.polygon = [h.polygon; [x y]];
    if size(h.polygon,1)>1
        line(h.polygon(end-1:end,1),h.polygon(end-1:end,2),'color','w')
    end
    guidata(hObject,h);
end

function update_display_mode(hObject, h)
switch h.display_mode.Value
    case 1 % mean image - green
        h.dat.map = 2;
        h.dat.procmap = 0;
        redraw_meanimg_PG(h);
    case 2 % mean image - red
        h.dat.map = 3;
        redraw_meanimg_PG(h);
    case 3 % mean image - proc
        h.dat.map = 2;
        %h.dat.procmap = 1 -  h.dat.procmap;
        h.dat.procmap = 1;
        if h.dat.map>1
            redraw_meanimg_PG(h);
        end
    case 4 % ROIs only
        h.dat.map = 1;
        redraw_figure(h);
    case 5 % overlay included ROIs
        h.dat.map = 2;
        redraw_meanimg_overlay(h);
    case 6 % overlay all ROIs
        h.dat.map = 2;
        redraw_meanimg_overlay(h,1);
    case 7 % overlay external ROIs
        % get ROI mask image from independent session
        if isempty(h.ROIs_compare)
            [filename1,filepath1]=uigetfile(fileparts(h.dat.filename), 'Select ROI image to compare');
            if filename1
                h.ROIs_compare = load(fullfile(filepath1, filename1));
            end
        end
        if ~isempty(h.ROIs_compare)
            guidata(hObject,h);
            h.dat.map = 2;
            redraw_meanimg_overlay(h,2);
        end
    case 8 % ROIs vs ROIs
end

guidata(hObject,h);

% --- Executes on selection change in display_mode.
function display_mode_Callback(hObject, eventdata, h)
% hObject    handle to display_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns display_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from display_mode
update_display_mode(hObject, h)

function mask_multiplier_1_Callback(hObject, eventdata, h)
% hObject    handle to mask_multiplier_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ismember(get(h.display_mode,'Value'),[5,6,7])
    update_display_mode(hObject, h)
end

function mask_multiplier_2_Callback(hObject, eventdata, h)
% hObject    handle to mask_multiplier_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ismember(get(h.display_mode,'Value'),[5,6,7])
    update_display_mode(hObject, h)
end

function row_offset_Callback(hObject, eventdata, h)
% hObject    handle to row_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(h.display_mode,'Value')==7) && ~isempty(h.ROIs_compare)
    update_display_mode(hObject, h)
end

function column_offset_Callback(hObject, eventdata, h)
% hObject    handle to column_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(h.display_mode,'Value')==7) && ~isempty(h.ROIs_compare)
    update_display_mode(hObject, h)
end

function update_F_display_Callback(hObject, eventdata, h)
[h] = update_cell_neuropil_traces(h);
redraw_fluorescence_PG(h)
guidata(hObject,h);


% --- Executes on button press in define_polygon.
function define_polygon_Callback(hObject, eventdata, h)
% hObject    handle to define_polygon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if h.define_polygon.Value
    h.define_polygon.BackgroundColor = 'red';
    h.polygon = [];
else
    h.define_polygon.BackgroundColor = 'black';
    h.polygon = [h.polygon; h.polygon(1,:)];
end
guidata(hObject,h);
% Hint: get(hObject,'Value') returns toggle state of define_polygon


% --- Executes on button press in show_polygon.
function show_polygon_Callback(hObject, eventdata, h)
% hObject    handle to show_polygon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(h.polygon)
    if h.show_polygon.Value
        h.show_polygon.BackgroundColor = 'red';
        polygon_overlay(h);
    else
        h.show_polygon.BackgroundColor = 'black';
        update_display_mode(hObject, h);
    end
end
% Hint: get(hObject,'Value') returns toggle state of show_polygon


% --- Executes on button press in exclude_polygon.
function exclude_polygon_Callback(hObject, eventdata, h)
% hObject    handle to exclude_polygon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(h.polygon) && h.show_polygon.Value
    size(h.dat.res.iclust)
    my_polygon = poly2mask(h.polygon(:,1), h.polygon(:,2), ...
        size(h.dat.res.iclust,1), size(h.dat.res.iclust,2));
    clust_indexes = unique(h.dat.res.iclust(find(my_polygon)));
    h.dat.cl.manual(clust_indexes) = -0.5;
    h = splitROIleftright_PG(h);
    Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
    Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
    h.dat.img1.Sat     = Sat;
    h.dat.img2.Sat     = Sat;
    h = buildLambdaValue(h);
    update_display_mode(hObject, h); %redraw_figure(h);
    guidata(hObject,h);
end
% Hint: get(hObject,'Value') returns toggle state of exclude_polygon


% --- Executes on button press in rescale_mean_image.
function rescale_mean_image_Callback(hObject, eventdata, h)
% hObject    handle to rescale_mean_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_display_mode(hObject, h); %redraw_figure(h);
guidata(hObject,h);
% Hint: get(hObject,'Value') returns toggle state of rescale_mean_image
