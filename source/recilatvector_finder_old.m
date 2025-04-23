function varargout = recilatvector_finder(varargin)
% RECILATVECTOR_FINDER MATLAB code for recilatvector_finder.fig
%      RECILATVECTOR_FINDER, by itself, creates a new RECILATVECTOR_FINDER or raises the existing
%      singleton*.
%
%      H = RECILATVECTOR_FINDER returns the handle to a new RECILATVECTOR_FINDER or the handle to
%      the existing singleton*.
%
%      RECILATVECTOR_FINDER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RECILATVECTOR_FINDER.M with the given input arguments.
%
%      RECILATVECTOR_FINDER('Property','Value',...) creates a new RECILATVECTOR_FINDER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before recilatvector_finder_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to recilatvector_finder_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help recilatvector_finder

% Last Modified by GUIDE v2.5 04-Jan-2025 23:31:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @recilatvector_finder_OpeningFcn, ...
                   'gui_OutputFcn',  @recilatvector_finder_OutputFcn, ...
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


% --- Executes just before recilatvector_finder is made visible.
function recilatvector_finder_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to recilatvector_finder (see VARARGIN)

% Choose default command line output for recilatvector_finder
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
if numel(varargin)>1
    if strcmp(get(varargin{1}, 'type'), 'uimenu')
        if strcmp(varargin{1}.Text, 'Reciprocal Lattice Vector Tool')
            pbsetdatafigure_Callback(hObject, gcbf, handles)
            ax = findobj(gcbf, 'type', 'axes');
            set(ax, 'Ydir', 'normal');
        end
    end
end
init_spacegroup(hObject, handles)
% UIWAIT makes recilatvector_finder wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function init_spacegroup(varargin)
[sg, sgN, syst] = spacegrouplist;
sg = ['P', sg];
sgN = [0, sgN];
syst = ['none', syst];
nsg = {};
sg_all = load('sg_all.mat');
%set(gcf, 'userdata', sg_all');
for i=1:numel(sg)
    nsg{i} = sprintf('%i, %s: %s', sgN(i), syst{i}, sg{i});
end
hobj = varargin{1};
handles = varargin{2};
set(hobj, 'userdata', sg_all);
set(handles.pm_centering_type, 'string', nsg);
set(handles.pm_centering_type, 'value', 1);



% --- Outputs from this function are returned to the command line.
function varargout = recilatvector_finder_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pbloadcellinfo.
function pbloadcellinfo_Callback(hObject, eventdata, handles)
% hObject    handle to pbloadcellinfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%ct = evalin('base', 'cellinfo');

ct = getappdata(gcbf, 'cellinfo');
if isempty(ct)
    fprintf('Cellinfo has not been loaded. Therefore, loading from workspace.\n');
    ct = evalin('base', 'cellinfo');
end
fh = getappdata(gcbf, 'dataFigure');
delete(findobj(fh, 'type', 'line'));

if handles.ch_isReciprocal.Value
    for i=1:3
        set(handles.(['rbview',num2str(i)]), 'Value', 1);
        pbsetview_Callback(ct.recilatticevectors(i, :), [], handles);
    end

    lh = drawunitcell(ct.recilatticevectors(1, :), ...
        ct.recilatticevectors(2, :), ...
        ct.recilatticevectors(3, :));
else
    for i=1:3
        set(handles.(['rbview',num2str(i)]), 'Value', 1);
        pbsetview_Callback(ct.latticevectors(i, :), [], handles);
    end

    lh = drawunitcell(ct.latticevectors(1, :), ...
        ct.latticevectors(2, :), ...
        ct.latticevectors(3, :));
end
for i=1:numel(lh)
    set(lh, 'linestyle', ':', 'tag', 'unitcell');
end
if handles.ch_isReciprocal.Value
    set(handles.edview1, 'string', norm(ct.recilatticevectors(1, :)));
    set(handles.edview2, 'string', norm(ct.recilatticevectors(2, :)));
    set(handles.edview3, 'string', norm(ct.recilatticevectors(3, :)));
else
    set(handles.edview1, 'string', norm(ct.latticevectors(1, :)));
    set(handles.edview2, 'string', norm(ct.latticevectors(2, :)));
    set(handles.edview3, 'string', norm(ct.latticevectors(3, :)));
end
% --- Executes on button press in pbsetview.
function pbsetview_Callback(hObject, eventdata, handles)
% hObject    handle to pbsetview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ishandle(hObject)
    as = get_view;
else
    as = hObject;
end
switch lower(handles.bgview.SelectedObject.String)
    case 'view2'
        set(handles.edview2, 'string', norm(as));
        setappdata(gcbf, 'bs', as);
%         setappdata(gcbf, 'bRmat', Rmat);
        c = 'g';
        linetag = 'lineb';
    case 'view1'
        set(handles.edview1, 'string', norm(as));
        setappdata(gcbf, 'as', as);
%         setappdata(gcbf, 'aRmat', Rmat);
        c = 'r';
        linetag = 'linea';
    case 'view3'
        set(handles.edview3, 'string', norm(as));
        setappdata(gcbf, 'cs', as);
%         setappdata(gcbf, 'cRmat', Rmat);
        c = 'b';
        linetag = 'linec';
end
drawline(as, linetag, c, handles);

function h = drawline(v, linetag, col, handles)
fh = getappdata(gcbf, 'dataFigure');
ax = findobj(fh, 'type', 'axes');
xl = ax.XLim;
yl = ax.YLim;
zl = ax.ZLim;
ud = get(gcbf, 'userdata');
if isfield(ud, 'cm')
    cm = ud.cm;
else
    cm = [0, 0, 0];
end
if handles.ch_isReciprocal.Value
    cm = [0,0,0];
else
    % if cm is out of view
    if cm(1)<xl(1) | cm(1)>xl(2) | cm(2)<yl(1) | cm(2)>yl(2) | cm(3)<zl(1) | cm(3)>zl(2)
        cm(1) = mean(xl);
        cm(2) = mean(yl);
        cm(3) = mean(zl);
        ud.cm = cm;
        set(gcbf, 'userdata', ud);
    end
end
if isempty(fh)
    error('dataFigure has not been selected yet.')
end
%ax = findobj(fh, 'type', 'axes');
delete(findobj(ax, 'tag', linetag));
h = line(ax, [0, v(1)]+cm(1), [0, v(2)]+cm(2), [0, v(3)]+cm(3));
set(h, 'color', col);
set(h, 'tag', linetag);

% --- Executes on button press in pbcleargrid.
function pbcleargrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbcleargrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fh = getappdata(gcbf, 'dataFigure');
ax = findobj(fh, 'type', 'axes');
delete(findobj(fh, 'tag', 'grid_hkl'));
delete(findobj(fh, 'tag', 'grid_hkl_P'));
set(findobj(ax, 'type', 'patch'), 'facealpha', 1);
if handles.mn_overlay2D.Checked
    fh = getappdata(gcbf, 'expsumimgFigure');
    if ishandle(fh)
        delete(findobj(fh, 'tag', 'simulated_hkls'));
    end
end


% --- Executes on button press in pbview3.
function pbview3_Callback(hObject, eventdata, handles)
% hObject    handle to pbview3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edview1_Callback(hObject, eventdata, handles)
% hObject    handle to edview1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edview1 as text
%        str2double(get(hObject,'String')) returns contents of edview1 as a double
set_recivector(hObject, handles)

function set_recivector(hObject, handles)
    val = str2double(get(hObject,'String'));
    switch lower(hObject.Tag)
        case 'edview1'
            vname = 'a';
        case 'edview2'
            vname = 'b';
        case 'edview3'
            vname = 'c';
    end
    if handles.ch_isReciprocal.Value
        vname = [vname, 's'];
    end
    v = getappdata(gcbf, vname);
    v = v/norm(v)*val;
    setappdata(gcbf, vname, v);
    fh = getappdata(gcbf, 'dataFigure');
    vh = findobj(fh, 'tag', ['line', vname(1)]);
    if ~isempty(vh)
        set_length(v, vh);
    end
    construct_cellinfo(handles)
    %update_from_cellinfo(handles)
    update_grid(handles)

function set_length(v, vhandle)
    udf = get(gcbf, 'userdata');
    if isfield(udf, 'cm')
        cm = udf.cm;
    else
        cm = [0,0,0];
    end
    if contains(vhandle.Tag, 's')
        cm = [0,0,0];
    end
    set(vhandle, 'xdata', [0, v(1)]+cm(1));
    set(vhandle, 'ydata', [0, v(2)]+cm(2));
    set(vhandle, 'zdata', [0, v(3)]+cm(3));

% --- Executes during object creation, after setting all properties.
function edview1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edview1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edview2_Callback(hObject, eventdata, handles)
% hObject    handle to edview2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edview2 as text
%        str2double(get(hObject,'String')) returns contents of edview2 as a double
set_recivector(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edview2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edview2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edview3_Callback(hObject, eventdata, handles)
% hObject    handle to edview3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edview3 as text
%        str2double(get(hObject,'String')) returns contents of edview3 as a double
set_recivector(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edview3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edview3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slview1_Callback(hObject, eventdata, handles)
% hObject    handle to slview1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_callback(handles);



% --- Executes during object creation, after setting all properties.
function slview1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slview1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slview2_Callback(hObject, eventdata, handles)
% hObject    handle to slview2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_callback(handles)

% --- Executes during object creation, after setting all properties.
function slview2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slview2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slview3_Callback(hObject, eventdata, handles)
% hObject    handle to slview3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_callback(handles)

% --- Executes during object creation, after setting all properties.
function slview3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slview3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pbprint.
function pbprint_Callback(hObject, eventdata, handles)
% hObject    handle to pbprint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    construct_cellinfo(handles)
    ct = getappdata(handles.figure1, 'cellinfo');

    if handles.ch_isReciprocal.Value
        fprintf('\n');
        fprintf('Reciprocal lattice vectors ==========\n')
        fprintf('a* = [%0.4e, %0.4e, %0.4e]\n', ct.recilatticevectors(1, :));
        fprintf('b* = [%0.4e, %0.4e, %0.4e]\n', ct.recilatticevectors(2, :));
        fprintf('c* = [%0.4e, %0.4e, %0.4e]\n', ct.recilatticevectors(3, :));
    
        fprintf('Lattice vectors ==========\n')
        fprintf('a = [%0.4e, %0.4e, %0.4e]\n', ct.latticevectors(1, :));
        fprintf('b = [%0.4e, %0.4e, %0.4e]\n', ct.latticevectors(2, :));
        fprintf('c = [%0.4e, %0.4e, %0.4e]\n', ct.latticevectors(3, :));
    
        fprintf('Lattice Parameters ==========\n')
        fprintf('[%0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f]\n', ...
            ct.A, ct.B, ct.C, ct.alpha, ct.beta, ct.gamma);
        fprintf('Volume = %0.4e\n', ct.Vol);
    else
        fprintf('\n');
        fprintf('Lattice vectors ==========\n')
        fprintf('a = [%0.4e, %0.4e, %0.4e]\n', ct.latticevectors(1, :));
        fprintf('b = [%0.4e, %0.4e, %0.4e]\n', ct.latticevectors(2, :));
        fprintf('c = [%0.4e, %0.4e, %0.4e]\n', ct.latticevectors(3, :));
    
        fprintf('Reciprocal lattice vectors ==========\n')
        fprintf('a* = [%0.4e, %0.4e, %0.4e]\n', ct.recilatticevectors(1, :));
        fprintf('b* = [%0.4e, %0.4e, %0.4e]\n', ct.recilatticevectors(2, :));
        fprintf('c* = [%0.4e, %0.4e, %0.4e]\n', ct.recilatticevectors(3, :));
    
        fprintf('Lattice Parameters ==========\n')
        fprintf('[%0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f]\n', ...
            ct.A, ct.B, ct.C, ct.alpha, ct.beta, ct.gamma);
        fprintf('Volume = %0.4e\n', ct.Vol);
    end
%    setappdata(gcbf, 'cellinfo', ct)
%    update_from_cellinfo(handles)
    assignin('base', 'cellinfo', ct);

function construct_cellinfo(handles)
    
    if handles.ch_isReciprocal.Value
        as = getappdata(gcbf, 'as');
        bs = getappdata(gcbf, 'bs');
        cs = getappdata(gcbf, 'cs');
        ct = celcon4rcplatvec(as, bs, cs);
    else
        a = getappdata(gcbf, 'a');
        b = getappdata(gcbf, 'b');
        c = getappdata(gcbf, 'c');
        ct = celcon4latvec(a, b, c);
    end
    setappdata(gcbf, 'cellinfo', ct)
    update_from_cellinfo(handles)

% --- Executes on button press in pbgrid.
function pbgrid_Callback(hObject, eventdata, handles)
% hObject    handle to pbgrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pbcleargrid_Callback([], [], handles);
ud = get(gcbf, 'userdata');
if isfield(ud, 'cm')
    cm = ud.cm;
else
    cm = [0, 0, 0];
end

fh = getappdata(gcbf, 'dataFigure');
if handles.ch_isReciprocal.Value
    v1 = getappdata(gcbf, 'as');
    v2 = getappdata(gcbf, 'bs');
    v3 = getappdata(gcbf, 'cs');
else
    v1 = getappdata(gcbf, 'a');
    v2 = getappdata(gcbf, 'b');
    v3 = getappdata(gcbf, 'c');
end
%ev1 = findobj(gcbf, 'tag', 'edview1');
%ev2 = findobj(gcbf, 'tag', 'edview2');
%ev3 = findobj(gcbf, 'tag', 'edview3');
%alen = str2double(get(ev1, 'String'));
%blen = str2double(get(ev2, 'String'));
%clen = str2double(get(ev3, 'String'));

ax = findobj(fh, 'type', 'axes');
xl = (get(ax, 'XLim'));
yl = (get(ax, 'YLim'));
zl = (get(ax, 'ZLim'));
box = [xl(1), yl(1), zl(1);xl(2),yl(1),zl(1);
    xl(1), yl(2), zl(1);xl(2),yl(2),zl(1);
    xl(1),yl(2),zl(2);xl(2),yl(2),zl(2);
    xl(1),yl(1),zl(2);xl(2),yl(1),zl(2);];
box = box-cm;
latticevectors = [v1;v2;v3];
% fractioanl coordinates
hbox = box/latticevectors;
hmin = round(min(hbox(:,1)));
hmax = round(max(hbox(:,1)));
kmin = round(min(hbox(:,2)));
kmax = round(max(hbox(:,2)));
lmin = round(min(hbox(:,3)));
lmax = round(max(hbox(:,3)));
pmcobj = findobj(gcbf, 'tag', 'pm_centering_type');
str = get(pmcobj, 'String');
strv = get(pmcobj, 'Value');
%hmax = max(abs(fix(xl/alen)))+2;
%kmax = max(abs(fix(yl/blen)))+2;
%lmax = max(abs(fix(zl/clen)))+2;
msobj = findobj(gcbf, 'tag', 'ed_markersize');
ms = str2double(get(msobj, 'String'));

hl = linspace(hmin, hmax, hmax-hmin+1);
kl = linspace(kmin, kmax, kmax-kmin+1);
ll = linspace(lmin, lmax, lmax-lmin+1);

[X, Y, Z] = meshgrid(hl, kl, ll);
hkls = [X(:), Y(:), Z(:)];
N = strtok(str{strv}, ',');
N = str2double(N);
switch N
    case 0
        hkls_cen = [];
        sg = [];
    otherwise
        if handles.ch_isReciprocal.Value
            ud = get(gcbf, 'userdata');
            sg = ud.sg{N};
            hkls_cen = hkl_listall(sg, [hmax, kmax, lmax]);
        else
            sg = 'I';
            hkls_cen = hkls + [0.5, 0.5, 0.5];
        end
end
% convert fractional to cartesian coordinate
np = hkls*latticevectors;
if ~handles.ch_isReciprocal.Value
    np = np+cm;
end
% remove peaks out of the view
t = np(:,1)<xl(1) | np(:,1)>xl(2);
t = t | np(:,2)<yl(1) | np(:,2)>yl(2);
t = t | np(:,3)<zl(1) | np(:,3)>zl(2);
np(t, :) = [];
hkls(t, :) = [];

qmax = str2double(handles.grid_q_max.String);
d = sqrt(sum(np.^2, 2));
t = d>qmax;
np(t, :) = [];
hkls(t, :) = [];

if ~isempty(hkls_cen)
    np_cen = hkls_cen*latticevectors;
    if ~handles.ch_isReciprocal.Value
        np_cen = np_cen+cm;
    end
    d = sqrt(sum(np_cen.^2, 2));
    t = d>qmax;
    hkls_cen(t, :) = [];
end
fh = getappdata(gcbf, 'dataFigure');
ud = get(fh, 'userdata');
if ~isempty(sg)
    ud.sginfo = sg;
end
if ~handles.ch_isReciprocal.Value
    if ~isfield(ud, 'Xd')
        ud.Xd = ud.X*ud.pixelsize;
        ud.Yd = ud.Y*ud.pixelsize;
        ud.Zd = ud.Z*ud.pixelsize;
    %    set(fh, 'userdata', ud);
    end
end
set(fh, 'userdata', ud);
set(findobj(ax, 'type', 'patch'), 'facealpha', 0.2);



if ~isempty(hkls_cen)
    h = plot3(ax, np(:,1), np(:,2), np(:,3), 'bo', 'markersize', 4, 'markerfacecolor', 'k');
    set(h, 'tag', 'grid_hkl_P');
    %set(h, 'userdata', hkls);
    cbh = findobj(gcbf, 'tag', 'cb_draw_P');
    if get(cbh, 'value')
        set(h, 'visible', 'on')
    else
        set(h, 'visible', 'off')
    end
    np = hkls_cen*latticevectors;
    t = np(:,1)<xl(1) | np(:,1)>xl(2);
    t = t | np(:,2)<yl(1) | np(:,2)>yl(2);
    t = t | np(:,3)<zl(1) | np(:,3)>zl(2);
    np(t, :) = [];
    if ~handles.ch_isReciprocal.Value
        np = np+cm;
    end
    hkls_cen(t, :) = [];
    %intv = interp3(ud.Xax, ud.Yax, ud.Zax,ud.map, np(:,1),np(:,2),np(:,3));
    %t = isnan(intv);
    h = plot3(ax, np(:,1), np(:,2), np(:,3), 'ro', 'markersize', ms);
    %h = plot3(ax, np(t,1), np(t,2), np(t,3), 'ro', 'markersize', ms);
    set(h, 'tag', 'grid_hkl');
    set(h, 'userdata', hkls_cen);
    %h = plot3(ax, np(~t,1), np(~t,2), np(~t,3), 'ro', 'markersize', ms, 'markerfacecolor', 'r');
    %set(h, 'tag', 'grid_hkl');
    %set(h, 'userdata', hkls_cen);
else
     if ~handles.ch_isReciprocal.Value
        try
            intv = interp3(ud.Xd, ud.Yd, ud.Zd,ud.map, np(:,1),np(:,2),np(:,3));
        catch
            intv = interp3(ud.Yd, ud.Xd, ud.Zd,ud.map, np(:,2),np(:,1),np(:,3));
        end
        t = isnan(intv);
     else
         t = 1:size(np, 1);
     end
    h = plot3(ax, np(t,1), np(t,2), np(t,3), 'ro', 'markersize', ms);
    set(h, 'tag', 'grid_hkl');
    set(h, 'userdata', hkls_cen);
    h = plot3(ax, np(~t,1), np(~t,2), np(~t,3), 'ro', 'markersize', ms, 'markerfacecolor', 'r');
    set(h, 'tag', 'grid_hkl');
    set(h, 'userdata', hkls_cen);
end


if handles.mn_overlay2D.Checked
    fh = getappdata(gcbf, 'expsumimgFigure');
    if ishandle(fh)
        calculate_q_pixel
    end
end


% --- Executes on button press in pblclearlines.
function pblclearlines_Callback(hObject, eventdata, handles)
% hObject    handle to pblclearlines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fh = getappdata(gcbf, 'dataFigure');
lh = findobj(fh, 'type', 'line', '-and', '-not', {'Tag', ''});
k = [];
for i=1:numel(lh)
    if contains(lh(i).Tag, 'line')
        k = [k, i];
    end
end
delete(lh(k));
    
% --- Executes on button press in pbgrid3.
function pbgrid3_Callback(hObject, eventdata, handles)
% hObject    handle to pbgrid3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function [as, Rmat] = get_view(varargin)
fh = getappdata(gcbf, 'dataFigure');
if isempty(fh)
    error('Choose a figure with a map of interest by pressing "Set Data Figure".')
end
[as, Rmat] = view2vect(fh);

% Decide the length of as based on the box size..
ax = findobj(fh, 'type', 'axes');
zl = get(ax, 'zlim');
yl = get(ax, 'ylim');
xl = get(ax, 'xlim');
maxl = max(sqrt(xl.^2+yl.^2+zl.^2));
as = as*maxl;

% Get the rotation matrix.
Rmat = Rmat(1:3, 1:3); % Rmat rotates v to [0, 0, 1]

%as = viewvector/norm(viewvector);


% --- Executes when selected object is changed in bgview.
function bgview_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in bgview 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fh = getappdata(gcbf, 'dataFigure');
switch lower(handles.bgview.SelectedObject.String)
    case 'view2'
        name = 'b';
    case 'view1'
        name = 'a';
    case 'view3'
        name = 'c';
end

if handles.ch_isReciprocal.Value
    name = [name, 's'];
end

v = getappdata(gcbf, name);

if ~isempty(v)
    figure(fh)
    view(v)
end


% --- Executes when selected object is changed in bgview.
function bgview_cross_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in bgview 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fh = getappdata(gcbf, 'dataFigure');
switch lower(handles.bgview_cross.SelectedObject.String)
    case 'v2 x v3'
        name = 'a';
    case 'v3 x v1'
        name = 'b';
    case 'v1 x v2'
        name = 'c';
end

if ~handles.ch_isReciprocal.Value
    name = [name, 's'];
end

v = getappdata(gcbf, name);

if ~isempty(v)
    figure(fh)
    view(v)
end

                


% --- Executes on button press in pbsetdatafigure.
function pbsetdatafigure_Callback(varargin)
% hObject    handle to pbsetdatafigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ishandle(varargin{2})
    setappdata(varargin{1}, 'dataFigure', varargin{2});
else
    setappdata(gcbf, 'dataFigure', gcf);
end
handles = varargin{3};
set(handles.uipanel5, 'visible', 'off');
if handles.ch_isReciprocal.Value
    xlabel('q_x')
    ylabel('q_y')
    zlabel('q_z')
else
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
end
hold on
%axis image;

function pbsetexpsumfigure_Callback(varargin)
% hObject    handle to pbsetdatafigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ishandle(varargin{2})
    setappdata(varargin{1}, 'expsumimgFigure', varargin{2});
else
    setappdata(gcbf, 'expsumimgFigure', gcf);
end
%handles = varargin{3};
%set(handles.uipanel5, 'visible', 'off');
xlabel('X')
ylabel('Y')
hold on
axis image;axis xy;


function edshiftcenterval_Callback(hObject, eventdata, handles)
% hObject    handle to edshiftcenterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edshiftcenterval as text
%        str2double(get(hObject,'String')) returns contents of edshiftcenterval as a double

% --- Executes during object creation, after setting all properties.
function edshiftcenterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edshiftcenterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbminus.
function pbminus_Callback(hObject, eventdata, handles)
% hObject    handle to pbminus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ds = str2double(handles.edshiftcenterval.String);
shiftcenter(handles, -1, ds);


% --- Executes on button press in pbplus.
function pbplus_Callback(hObject, eventdata, handles)
% hObject    handle to pbplus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ds = str2double(handles.edshiftcenterval.String);
shiftcenter(handles, 1, ds);

function shiftcenter(handles, sign, ds)
switch lower(handles.bgshiftcenter.SelectedObject.String)
    case 'axis 1'
        vname = 'a';
    case 'axis 2'
        vname = 'b';
    case 'axis 3'
        vname = 'c';
end
if handles.ch_isReciprocal.Value
    vname = [vname, 's'];
end
v = getappdata(gcbf, vname);
fh = getappdata(gcbf, 'dataFigure');
%ax = findobj(fh, 'type', 'axes');
%lh = findobj(ax, 'tag', 'grid_hkl');
%lh = findobj(ax, 'type', 'line');
ud = get(gcbf, 'userdata');
if sign < 0
    ds = ds*(-1);
end
if ~isempty(ud)

    if handles.ch_isReciprocal.Value
        obj = findobj(fh, 'type', 'patch');
        obj.Vertices = obj.Vertices - v*ds;
        obj.XData = obj.XData - v(1)*ds;
        obj.YData = obj.YData - v(2)*ds;
        obj.ZData = obj.ZData - v(3)*ds;
    else
        if isfield(ud, 'cm')
            cm = ud.cm;
        else
            cm = [0,0,0];
        end
        cm = cm+v*ds;
        ud.cm = cm;
        fprintf('cm is at [%0.3f, %0.3f, %0.3f]\n', cm)
        set(gcbf, 'userdata', ud)
        pbdrawlines_Callback([],[],handles);
    end
end

% --- Executes on button press in pbdrawlines.
function pbdrawlines_Callback(hObject, eventdata, handles)
% hObject    handle to pbdrawlines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
draw_latvec(1, handles);
draw_latvec(2, handles);
draw_latvec(3, handles);

function draw_latvec(n, handles)
switch n
    case {1, 'as', 'a'}
        name = 'a';
        c = 'r';
        linetag = 'linea';
    case {2, 'bs', 'b'}
        name = 'b';
        c = 'g';
        linetag = 'lineb';
    case {3, 'cs', 'c'}
        name = 'c';
        c = 'b';
        linetag = 'linec';
end
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
v = getappdata(handles.figure1, name);
drawline(v, linetag, c, handles);
update_grid(handles)

function update_grid(varargin)
    handles = varargin{1};
    fh = getappdata(gcbf, 'dataFigure');
    h = findobj(fh, 'tag', 'grid_hkl');
    if ~isempty(h)
        delete(h);
        pbgrid_Callback([], [], handles);
    end

% --- Executes on button press in pbobservedhkl.
function pbobservedhkl_Callback(hObject, eventdata, handles)
% hObject    handle to pbobservedhkl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
as = getappdata(gcbf, 'as');
bs = getappdata(gcbf, 'bs');
cs = getappdata(gcbf, 'cs');

fh = getappdata(gcbf, 'dataFigure');
ax = findobj(fh, 'type', 'axes');
%lh = findobj(ax, 'type', 'line', '-and', {'Tag', ''});
lh = findobj(ax, 'type', 'line');
if isempty(lh)
    fprintf('Data Figure does not have any "line" object.\n')
    return
end
for i=1:numel(lh)
    l = lh(i);
    xd = get(l, 'xdata');
    if numel(xd) == 1
        continue
    end
    yd = get(l, 'ydata');
    zd = get(l, 'zdata');
end    
vv = [xd(:), yd(:), zd(:)];

ohkl = check_hit_hkl(as, bs, cs, vv);
fprintf('\n');
fprintf('List of observed hkls ==========\n');
for i=1:size(ohkl, 1)
    fprintf('[%i, %i, %i]\n', ohkl(i,:))
end
fprintf('Total %i reflections are observed.\n', size(ohkl, 1));


function slider_callback(varargin)
handles = varargin{1};
sv = get(gcbo, 'value');
sv = sv - 0.5;
switch get(gcbo, 'Tag')
    case 'slview1'
        name = 'a';
    case 'slview2'
        name = 'b';
    case 'slview3'
        name = 'c';
end
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
v = getappdata(gcbf, name);
vp = get_view;
v = rotate_around_vector(v, vp, sv*90);
setappdata(gcbf, name, v);

construct_cellinfo(handles)

draw_latvec(name, handles);
set(gcbo, 'value', 0.5);


function slider_len_callback(handles)
sv = get(gcbo, 'value');
sv = sv - 0.5;
switch get(gcbo, 'Tag')
    case 'slalen'
        %name = 'as';
        n = str2double(handles.edview1.String);
        n = n*(1+sv);
        set(handles.edview1, 'string', n);
        edview1_Callback(handles.edview1, [], handles)
    case 'slblen'
        %name = 'bs';
        n = str2double(handles.edview2.String);
        n = n*(1+sv);
        set(handles.edview2, 'string', n);
        edview2_Callback(handles.edview2, [], handles)
    case 'slclen'
        %name = 'cs';
        n = str2double(handles.edview3.String);
        n = n*(1+sv);
        set(handles.edview3, 'string', n);
        edview3_Callback(handles.edview3, [], handles)
end

set(gcbo, 'value', 0.5);

% --- Executes on button press in pbseta.
function pbseta_Callback(hObject, eventdata, handles)
% hObject    handle to pbseta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cv = get_crossvector(handles);
name = 'a';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
a = getappdata(gcbf, name);
cv = cv/norm(cv);
if ~isempty(a)
    cv = cv*norm(a);
else
    set(handles.edview1, 'string', 1);
end

setappdata(gcbf, name, cv);
draw_latvec(1);

function cv = get_crossvector(handles)
v = get_view;
val = get(handles.lb_color,'Value');
switch val
    case 1
        name = 'a';
    case 2
        name = 'b';
    case 3
        name = 'c';
end
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
vt = getappdata(gcbf, name);
cv = cross(v, vt);

% --- Executes on button press in pbsetb.
function pbsetb_Callback(hObject, eventdata, handles)
% hObject    handle to pbsetb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cv = get_crossvector(handles);
name = 'bs';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
a = getappdata(gcbf, name);
cv = cv/norm(cv);
if ~isempty(a)
    cv = cv*norm(a);
else
    set(handles.edview2, 'string', 1);
end

setappdata(gcbf, name, cv);
draw_latvec(2);


% --- Executes on button press in pbsetc.
function pbsetc_Callback(hObject, eventdata, handles)
% hObject    handle to pbsetc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cv = get_crossvector(handles);
name = 'cs';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
a = getappdata(gcbf, name);
cv = cv/norm(cv);
if ~isempty(a)
    cv = cv*norm(a);
else
    set(handles.edview3, 'string', 1);
end
setappdata(gcbf, name, cv);
draw_latvec(3, handles);


% --- Executes on selection change in lb_color.
function lb_color_Callback(hObject, eventdata, handles)
% hObject    handle to lb_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_color


% --- Executes during object creation, after setting all properties.
function lb_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edsetview_Callback(hObject, eventdata, handles)
% hObject    handle to edsetview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edsetview as text
%        str2double(get(hObject,'String')) returns contents of edsetview as a double
vv = eval(get(hObject,'String'));
cellinfo = getappdata(gcbf, 'cellinfo');
if handles.ch_isReciprocal.Value
    vv = vv*cellinfo.recilatticevectors;
else
    vv = vv*cellinfo.latticevectors;
end
view(vv);
%pbsetview_Callback(vv, [], handles)

% --- Executes during object creation, after setting all properties.
function edsetview_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edsetview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbsetCc.
function pbsetCc_Callback(hObject, eventdata, handles)
% hObject    handle to pbsetCc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = getcursor;
set(handles.edview3, 'string', norm(a));
name = 'c';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
setappdata(gcbf, name, a);
draw_latvec(3);


% --- Executes on button press in pbsetCb.
function pbsetCb_Callback(hObject, eventdata, handles)
% hObject    handle to pbsetCb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = getcursor;
set(handles.edview2, 'string', norm(a));
name = 'b';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
setappdata(gcbf, name, a);
draw_latvec(2);


% --- Executes on button press in pbsetCa.
function pbsetCa_Callback(hObject, eventdata, handles)
% hObject    handle to pbsetCa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = getcursor;
set(handles.edview1, 'string', norm(a));
name = 'a';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
setappdata(gcbf, name, a);
draw_latvec(1);

function st = getcursor(varargin)
    dcm = datacursormode(gcf);
    s = getCursorInfo(dcm);
    if numel(s) > 1
        error('Choose one point')
    end
    st = s.Position;



function ed_markersize_Callback(hObject, eventdata, handles)
% hObject    handle to ed_markersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_markersize as text
%        str2double(get(hObject,'String')) returns contents of ed_markersize as a double
pbgrid_Callback([], [], handles)

% --- Executes during object creation, after setting all properties.
function ed_markersize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_markersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function m_File_Callback(hObject, eventdata, handles)
% hObject    handle to m_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function m_loadcellinfo_Callback(hObject, eventdata, handles)
% hObject    handle to m_loadcellinfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, pname] = uigetfile('*.mat', 'Pick a Map file');
if isequal(fname,0) || isequal(pname,0)
   disp('User pressed cancel')
else
   disp(['User selected ', fullfile(pname, fname)])
end
loadedData = load(fullfile(pname, fname));
if isfield(loadedData, 'cellinfo')
    cellinfo = loadedData.cellinfo;
else
    return
end
load_cellinfo(handles, cellinfo);
update_from_cellinfo(handles)


%assignin('base', 'cellinfo', cellinfo)
function load_cellinfo(handles, cellinfo)
    ev1 = findobj(handles.figure1, 'tag', 'edview1');
    ev2 = findobj(handles.figure1, 'tag', 'edview2');
    ev3 = findobj(handles.figure1, 'tag', 'edview3');
    if handles.ch_isReciprocal.Value
        as = cellinfo.recilatticevectors(1, :);
        bs = cellinfo.recilatticevectors(2, :);
        cs = cellinfo.recilatticevectors(3, :);
        aslen = norm(as);
        bslen = norm(bs);
        cslen = norm(cs);
    else
        a = cellinfo.latticevectors(1, :);
        b = cellinfo.latticevectors(2, :);
        c = cellinfo.latticevectors(3, :);
        aslen = norm(a);
        bslen = norm(b);
        cslen = norm(c);
    end

    if handles.ch_isReciprocal.Value
        setappdata(handles.figure1, 'as', as);
        setappdata(handles.figure1, 'bs', bs);
        setappdata(handles.figure1, 'cs', cs);
    else
        setappdata(handles.figure1, 'a', a);
        setappdata(handles.figure1, 'b', b);
        setappdata(handles.figure1, 'c', c);
    end
    if ~isfield(cellinfo, 'U')
            % calculate orientation matrix
        cp = [cellinfo.A,cellinfo.B,cellinfo.C,cellinfo.alpha,cellinfo.beta,cellinfo.gamma];
        t = celcon(cp);Bmat=t.recimat;
        UB = cellinfo.recimat;
        U = UB*inv(Bmat);
        cellinfo.U = U;
    end

    
    set(ev1, 'String', sprintf('%0.5f', aslen));
    set(ev2, 'String', sprintf('%0.5f', bslen));
    set(ev3, 'String', sprintf('%0.5f', cslen));
    
    setappdata(handles.figure1, 'cellinfo', cellinfo);

function update_from_cellinfo(handles)
    cellinfo = getappdata(handles.figure1, 'cellinfo');
    set(handles.ed_a, 'string', num2str(cellinfo.A))
    set(handles.ed_b, 'string', num2str(cellinfo.B))
    set(handles.ed_c, 'string', num2str(cellinfo.C))
    set(handles.ed_alpha, 'string', num2str(cellinfo.alpha))
    set(handles.ed_beta, 'string', num2str(cellinfo.beta))
    set(handles.ed_gamma, 'string', num2str(cellinfo.gamma))
    if isfield(cellinfo, 'U')
        [PHI, THETA, PSI] = rotation3dToEulerAngles(cellinfo.U);
    else
        cellinfo.U = eye(3);
        PHI = 0; THETA = 0; PSI = 0;
        setappdata(handles.figure1, 'cellinfo', cellinfo);
    end
    set(handles.ed_phi_euler, 'string', num2str(PHI))
    set(handles.ed_theta_euler, 'string', num2str(THETA))
    set(handles.ed_psi_euler, 'string', num2str(PSI))
    setappdata(handles.figure1, 'as', cellinfo.recilatticevectors(1, :));
    setappdata(handles.figure1, 'bs', cellinfo.recilatticevectors(2, :));
    setappdata(handles.figure1, 'cs', cellinfo.recilatticevectors(3, :));
    setappdata(handles.figure1, 'a', cellinfo.latticevectors(1, :));
    setappdata(handles.figure1, 'b', cellinfo.latticevectors(2, :));
    setappdata(handles.figure1, 'c', cellinfo.latticevectors(3, :));

    if handles.ch_isReciprocal.Value
        set(handles.edview1, 'String', sprintf('%0.5f', norm(cellinfo.recilatticevectors(1, :))));
        set(handles.edview2, 'String', sprintf('%0.5f', norm(cellinfo.recilatticevectors(2, :))));
        set(handles.edview3, 'String', sprintf('%0.5f', norm(cellinfo.recilatticevectors(3, :))));
    else
        set(handles.edview1, 'String', sprintf('%0.5f', norm(cellinfo.latticevectors(1, :))));
        set(handles.edview2, 'String', sprintf('%0.5f', norm(cellinfo.latticevectors(2, :))));
        set(handles.edview3, 'String', sprintf('%0.5f', norm(cellinfo.latticevectors(3, :))));
    end

function cellinfo = update_to_cellinfo(handles)
    A = str2double(handles.ed_a.String);
    B = str2double(handles.ed_b.String);
    C = str2double(handles.ed_c.String);
    alpha = str2double(handles.ed_alpha.String);
    beta = str2double(handles.ed_beta.String);
    gamma = str2double(handles.ed_gamma.String);
    phi = str2double(handles.ed_phi_euler.String);
    theta = str2double(handles.ed_theta_euler.String);
    psi = str2double(handles.ed_psi_euler.String);
    U = eulerAnglesToRotation3d(phi, theta, psi);
    cellinfo = celcon([A, B, C, alpha, beta, gamma], U(1:3, 1:3));
    setappdata(handles.figure1, 'cellinfo', cellinfo);
    setappdata(handles.figure1, 'as', cellinfo.recilatticevectors(1, :));
    setappdata(handles.figure1, 'bs', cellinfo.recilatticevectors(2, :));
    setappdata(handles.figure1, 'cs', cellinfo.recilatticevectors(3, :));
    setappdata(handles.figure1, 'a', cellinfo.latticevectors(1, :));
    setappdata(handles.figure1, 'b', cellinfo.latticevectors(2, :));
    setappdata(handles.figure1, 'c', cellinfo.latticevectors(3, :));
%     set(handles.edview1, 'String', sprintf('%0.5f', norm(cellinfo.recilatticevectors(1, :))));
%     set(handles.edview2, 'String', sprintf('%0.5f', norm(cellinfo.recilatticevectors(2, :))));
%     set(handles.edview3, 'String', sprintf('%0.5f', norm(cellinfo.recilatticevectors(3, :))));
    if handles.ch_isReciprocal.Value
        set(handles.edview1, 'String', sprintf('%0.5f', norm(cellinfo.recilatticevectors(1, :))));
        set(handles.edview2, 'String', sprintf('%0.5f', norm(cellinfo.recilatticevectors(2, :))));
        set(handles.edview3, 'String', sprintf('%0.5f', norm(cellinfo.recilatticevectors(3, :))));
    else
        set(handles.edview1, 'String', sprintf('%0.5f', norm(cellinfo.latticevectors(1, :))));
        set(handles.edview2, 'String', sprintf('%0.5f', norm(cellinfo.latticevectors(2, :))));
        set(handles.edview3, 'String', sprintf('%0.5f', norm(cellinfo.latticevectors(3, :))));
    end
function m_setDataFigure_Callback(hObject, eventdata, handles)
pbsetdatafigure_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function m_save_cellinfo_Callback(hObject, eventdata, handles)
% hObject    handle to m_save_cellinfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.ch_isReciprocal.Value
    as = getappdata(gcbf, 'as');
    bs = getappdata(gcbf, 'bs');
    cs = getappdata(gcbf, 'cs');
    cellinfo = celcon4rcplatvec(as, bs, cs);
else
    a = getappdata(gcbf, 'a');
    b = getappdata(gcbf, 'b');
    c = getappdata(gcbf, 'c');
    cellinfo = celcon4latvec(a, b, c);
end

[fname, pname] = uiputfile('*.mat', 'Save cellinfo As');
if isequal(fname,0) || isequal(pname,0)
   disp('User pressed cancel')
   return
end
fn = fullfile(pname, fname);
save(fn, 'cellinfo')
fprintf('Saved as %s.\n', fn)


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ishandle(hObject)
    as = get_view;
else
    as = hObject;
end

set(handles.edview1, 'string', norm(as));
name = 'a';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
setappdata(gcbf, name, as);
%         setappdata(gcbf, 'aRmat', Rmat);
c = 'r';
linetag = 'linea';

drawline(as, linetag, c, handles);

% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ishandle(hObject)
    as = get_view;
else
    as = hObject;
end
set(handles.edview2, 'string', norm(as));
name = 'b';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end

setappdata(gcbf, name, as);
%         setappdata(gcbf, 'bRmat', Rmat);
c = 'g';
linetag = 'lineb';

drawline(as, linetag, c, handles);

% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ishandle(hObject)
    as = get_view;
else
    as = hObject;
end

set(handles.edview3, 'string', norm(as));
name = 'c';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
setappdata(gcbf, name, as);
%         setappdata(gcbf, 'cRmat', Rmat);
c = 'b';
linetag = 'linec';

drawline(as, linetag, c, handles);


% --------------------------------------------------------------------
function m_draw_Callback(hObject, eventdata, handles)
% hObject    handle to m_draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function m_draw_unitcell_Callback(hObject, eventdata, handles)
% hObject    handle to m_draw_unitcell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ct = getappdata(gcbf, 'cellinfo');
if isempty(ct)
    fprintf('Cellinfo has not been loaded. Therefore, loading from workspace.\n');
    ct = evalin('base', 'cellinfo');
end
fh = getappdata(gcbf, 'dataFigure');
delete(findobj(fh, 'type', 'line'));
ax = findobj(fh, 'type', 'axes');
xl = get(ax, 'XLim');
yl = get(ax, 'YLim');
zl = get(ax, 'ZLim');
ud = get(gcbf, 'userdata');
if isfield(ud, 'cm')
    cm = ud.cm;
else
    cm = [0, 0, 0];
end
if handles.ch_isReciprocal.Value
    cm = [0,0,0];
end

fprintf('Origin of the unit cell is at [%0.4f, %0.4f, %0.4f]\n', cm);
vx = zeros(3,3);
for i=1:3
    set(handles.(['rbview',num2str(i)]), 'Value', 1);
    if handles.ch_isReciprocal.Value
        v = ct.recilatticevectors(i, :);
    else
        v = ct.latticevectors(i, :);
    end
    vx(i, 1:3) = v;
    pbsetview_Callback(v, [], handles);
end

lh = drawunitcell(vx(1, :), ...
    vx(2,:), ...
    vx(3,:));
for i=1:numel(lh)
    lh(i).XData = lh(i).XData+cm(1);
    lh(i).YData = lh(i).YData+cm(2);
    lh(i).ZData = lh(i).ZData+cm(3);
    set(lh, 'linestyle', ':', 'tag', 'unitcell');
end
set(ax, 'XLim', xl)
set(ax, 'YLim', yl)
set(ax, 'ZLim', zl)
% --------------------------------------------------------------------
function m_draw_EWsphere_Callback(hObject, eventdata, handles)
% hObject    handle to m_draw_EWsphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figh = getappdata(gcbf, 'dataFigure');
waveln = getappdata(figh, 'wavelength');
if isempty(waveln)
    eng = input('Input X-ray Energy (keV): ');
    waveln = eng2wl(eng);
    setappdata(figh, 'wavelength', waveln)
    fprintf('Type setappdata(gcf, "wavelength", [new value]) to replace.\n')
end
% list all posiible hkls
v1 = getappdata(gcbf, 'as');
v2 = getappdata(gcbf, 'bs');
v3 = getappdata(gcbf, 'cs');
ax = findobj(figh, 'type', 'axes');
xl = (get(ax, 'XLim'));
yl = (get(ax, 'YLim'));
zl = (get(ax, 'ZLim'));
[X, Y, Z] = meshgrid(-10:10, -10:10, -10:10);
np = [X(:), Y(:), Z(:)]*[v1; v2; v3];
t = np(:,1)<xl(1) | np(:,1)>xl(2);
t = t | np(:,2)<yl(1) | np(:,2)>yl(2);
t = t | np(:,3)<zl(1) | np(:,3)>zl(2);
np(t, :) = [];

% draw Ewald sphere.
cc = Ewaldsphere(waveln, np, 'y');
xd = get(cc, 'xdata');
yd = get(cc, 'ydata');
zd = get(cc, 'zdata');
dt{1} = xd;
dt{2} = yd;
dt{3} = zd;
set(cc, 'userdata', dt);
set(cc, 'facealpha', 0.9);
set(cc, 'tag', 'Ewaldsphere_on_RST');
set(handles.uipanel5, 'visible', 'on');
% --------------------------------------------------------------------
function m_clearall_Callback(hObject, eventdata, handles)
% hObject    handle to m_clearall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figh = getappdata(gcbf, 'dataFigure');
delete(findobj(figh, 'tag', 'Ewaldsphere_on_RST'));
delete(findobj(figh, 'tag', 'unitcell'));
set(handles.uipanel5, 'visible', 'off');

% --- Executes on slider movement.
function slider_phi_Callback(hObject, eventdata, handles)
% hObject    handle to slider_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
phi = get(hObject,'value');
set(handles.ed_phi_euler, 'string', num2str(phi));
theta = get(handles.slider_theta, 'value');
slider_rotate(phi, theta);

% --- Executes during object creation, after setting all properties.
function slider_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_theta_Callback(hObject, eventdata, handles)
% hObject    handle to slider_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
theta = get(hObject,'value');
set(handles.ed_theta_euler, 'string', theta);
phi = get(handles.slider_phi, 'value');
slider_rotate(phi, theta);

% --- Executes during object creation, after setting all properties.
function slider_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider_rotate(phi, theta)
figh = getappdata(gcbf, 'dataFigure');
cc = findobj(figh, 'tag', 'Ewaldsphere_on_RST');
if isempty(cc)
    fprintf('This function is to rotate Ewaldsphere. Draw the sphere to proceed.\n')
    return
end
dt = get(cc, 'userdata');
xd = dt{1};
yd = dt{2};
zd = dt{3};
R = eulerAnglesToRotation3d(phi, theta, 0);
R  = R(1:3, 1:3);
dt = [xd(:), yd(:), zd(:)]*R';
set(cc, 'xdata', reshape(dt(:,1), size(xd)));
set(cc, 'ydata', reshape(dt(:,2), size(xd)));
set(cc, 'zdata', reshape(dt(:,3), size(xd)));


function ed_phi_Callback(hObject, eventdata, handles)
% hObject    handle to ed_phi_euler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_phi_euler as text
%        str2double(get(hObject,'String')) returns contents of ed_phi_euler as a double
val = str2double(get(hObject, 'string'));
set(handles.slider_phi, 'value', val);
theta = str2double(get(handles.ed_theta_euler, 'string'));
slider_rotate(val, theta)

function ed_phi_euler_Callback(hObject, eventdata, handles)
% hObject    handle to ed_phi_euler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_phi_euler as text
%        str2double(get(hObject,'String')) returns contents of ed_phi_euler as a double
update_to_cellinfo(handles)


% --- Executes during object creation, after setting all properties.
function ed_phi_euler_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_phi_euler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function ed_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_phi_euler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function ed_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_phi_euler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_theta_Callback(hObject, eventdata, handles)
% hObject    handle to ed_theta_euler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_theta_euler as text
%        str2double(get(hObject,'String')) returns contents of ed_theta_euler as a double
val = str2double(get(hObject, 'string'));
set(handles.slider_theta, 'value', val);
phi = str2double(get(handles.ed_phi_euler, 'string'));
slider_rotate(phi, val)

function ed_theta_euler_Callback(hObject, eventdata, handles)
% hObject    handle to ed_theta_euler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_theta_euler as text
%        str2double(get(hObject,'String')) returns contents of ed_theta_euler as a double
update_to_cellinfo(handles)

% --- Executes during object creation, after setting all properties.
function ed_theta_euler_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_theta_euler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pm_centering_type.
function pm_centering_type_Callback(hObject, eventdata, handles)
% hObject    handle to pm_centering_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_centering_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_centering_type
pbgrid_Callback([], [], handles)

% --- Executes during object creation, after setting all properties.
function pm_centering_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_centering_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_draw_P.
function cb_draw_P_Callback(hObject, eventdata, handles)
% hObject    handle to cb_draw_P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_draw_P
figh = getappdata(gcbf, 'dataFigure');
cc = findobj(figh, 'tag', 'grid_hkl_P');
if isempty(cc)
    fprintf('Primitive grid was not plotted.\n')
    return
end

if get(hObject,'Value')
    set(cc, 'visible', 'on');
else
    set(cc, 'visible', 'off');
end


% --- Executes on button press in cb_selhkls.
function cb_selhkls_Callback(hObject, eventdata, handles)
% hObject    handle to cb_selhkls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_selhkls
ff = evalin('base', 'QmapFamily');
if isempty(ff)
    fprintf('No splitted patches yet.\n')
    return
end

fh = getappdata(gcbf, 'dataFigure');
ax = findobj(fh, 'type', 'axes');
ff = findobj(ax, 'type', 'patch');
if ~get(hObject,'Value')
    xl = (get(ax, 'XLim'));
    m = findobj(ax, 'tag', 'grid_hkl');
    set(ax, 'Units', 'Points')
    axpos = get(ax, 'Position');
    scale=diff(xl)/axpos(3)*m.MarkerSize;
    xd = m.XData;
    yd = m.YData;
    zd = m.ZData;
    pos = [xd(:), yd(:), zd(:)];
    for i=1:numel(ff)
        cm = mean(ff(i).Vertices);
        d = sqrt(sum((pos - cm).^2, 2));
        if any(d<scale)
            set(ff(i), 'visible', 'off')
        end
    end
else
    set(ff, 'visible', 'on');
end


% --- Executes on slider movement.
function slblen_Callback(hObject, eventdata, handles)
% hObject    handle to slblen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_len_callback(handles)

% --- Executes during object creation, after setting all properties.
function slblen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slblen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slclen_Callback(hObject, eventdata, handles)
% hObject    handle to slclen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_len_callback(handles)

% --- Executes during object creation, after setting all properties.
function slclen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slclen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slalen_Callback(hObject, eventdata, handles)
% hObject    handle to slalen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_len_callback(handles)

% --- Executes during object creation, after setting all properties.
function slalen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slalen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function m_edit_Callback(hObject, eventdata, handles)
% hObject    handle to m_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function m_autoindexing_Callback(hObject, eventdata, handles)
ff = evalin('base', 'QmapFamily');
if isempty(ff)
    warndlg('On the menu of draw_3dmap, run Tools>Split Isosurfaces.', 'Warning', 'modal');
    return
end

%fh = getappdata(gcbf, 'dataFigure');
%ax = findobj(fh, 'type', 'axes');
%ff = findobj(ax, 'type', 'patch');
fcm = zeros(numel(ff), 3);
for i=1:numel(ff)
    fcm(i, :) = ff{i}.cm;
end
ds = sqrt(sum(fcm.^2,2));
[~, ind1] = sort(ds);
if numel(ind1)>70
    fcm = fcm(ind1(1:70),:);
end
[lattice, indices] = index_dirax(fcm)

% 
%     set(handles.edview1, 'string', norm(a1));
%     setappdata(gcbf, 'as', a1);
%     c = 'r';
%     linetag = 'linea';
%     drawline(a1, linetag, c, handles);
% 
%     set(handles.edview2, 'string', norm(b1));
%     setappdata(gcbf, 'bs', b1);
%     c = 'g';
%     linetag = 'lineb';
%     drawline(b1, linetag, c, handles);
% 
%     set(handles.edview3, 'string', norm(c1));
%     setappdata(gcbf, 'cs', c1);
%     c = 'b';
%     linetag = 'linec';
%     drawline(c1, linetag, c, handles);
% % construct a unit cell.
%     pbprint_Callback([], [], handles);

function indexfound = find_allow_reflections(varargin)
handles = varargin{1};
ff = evalin('base', 'QmapFamily');
if isempty(ff)
    warndlg('On the menu of draw_3dmap, run Tools>Split Isosurfaces.', 'Warning', 'modal');
    return
end

fcm = zeros(numel(ff), 3);
for i=1:numel(ff)
    fcm(i, :) = ff{i}.cm;
end

fh = getappdata(gcbf, 'dataFigure');
if handles.ch_isReciprocal.Value
    v1 = getappdata(gcbf, 'as');
    v2 = getappdata(gcbf, 'bs');
    v3 = getappdata(gcbf, 'cs');
else
    v1 = getappdata(gcbf, 'a');
    v2 = getappdata(gcbf, 'b');
    v3 = getappdata(gcbf, 'c');
end

latticevectors = [v1;v2;v3];
al = norm(v1);
bl = norm(v2);
cl = norm(v3);
rmin = min([al,bl,cl])*0.1;
hkl = get(findobj(fh, 'tag', 'grid_hkl'), 'userdata');
qvs = hkl*latticevectors;
indexfound = [];
for i = 1:size(fcm,1)
    fv = fcm(i,:);
    ds = sqrt(sum((qvs - fv).^2,2));
    [dsmin, dsindx] = min(ds);
    if dsmin < rmin 
        hkl(dsindx, :)
        indexfound = [indexfound;hkl(dsindx,:)];
    end
end
assignin('base', 'indexfound', indexfound)


function m_autoindexing_Callback_old(hObject, eventdata, handles)
ff = evalin('base', 'QmapFamily');
if isempty(ff)
    warndlg('On the menu of draw_3dmap, run Tools>Split Isosurfaces.', 'Warning', 'modal');
    return
end

fh = getappdata(gcbf, 'dataFigure');
ax = findobj(fh, 'type', 'axes');
ff = findobj(ax, 'type', 'patch');
fcm = zeros(numel(ff), 3);
for i=1:numel(ff)
    fcm(i, :) = mean(ff(i).Vertices);
end

pmcobj = findobj(gcbf, 'tag', 'pm_centering_type');
str = get(pmcobj, 'String');
strv = get(pmcobj, 'Value');
hmax = 2;
kmax = 2;
lmax = 2;

% hl = linspace(-hmax, hmax, hmax*2+1);
% kl = linspace(-kmax, kmax, kmax*2+1);
% ll = linspace(-lmax, lmax, lmax*2+1);
% [X, Y, Z] = meshgrid(hl, kl, ll);
% hkls = [X(:), Y(:), Z(:)];
N = strtok(str{strv}, ',');
%sg = strtrim(sg(2:end));
N = str2double(N);
switch N
    case 0
        hkls_cen = [];
%     case 'I'
%         hkls_cen = [1, 1, 0; 1 -1, 0; 1, 0, 1];
%     case 'F'
%         hkls_cen = [1, 1, 1; 1, -1, -1; 1, -1, 1];
    otherwise
        ud = get(gcbf, 'userdata');
        hkls_cen = hkl_listall(ud.sg{N}, [hmax, kmax, lmax]);
end
if isempty(hkls_cen)
    refa = [1, 0, 0];
    refb = [0, 1, 0];
    refc = [0, 0, 1];
else
    refa = hkls_cen(1, :);
    k = 2;
    refb = hkls_cen(k, :);
%     while dot(refb, refa) ~= 0
%         k = k + 1;
%         refb = hkls_cen(k, :);
%     end
    refc = hkls_cen(k+1, :);
%     while (dot(refc, refa)~=0) | (dot(refc, refa)~=0)
%         k = k + 1;
%         refc = hkls_cen(k, :);
%     end
    m = [refa(:), refb(:), refc(:)];
    if det(m) == 0
        fprintf('Try different reference hkls so that DET(m) is not 0.\n')
    end
end

av = [1, 0, 0];
bv = [0, 1, 0];
cv = [0, 0, 1];
ang_range = [inf, inf, inf];


dis = sqrt(sum(fcm.^2, 2));

% determine a
    [dis, ind] = sort(dis, 'ascend');
    fcm = fcm(ind, :);
    % finding positive direction
    i = 1;
    a1 = fcm(i, :);

    phaseangle = angle2vect(av, fcm(1, :))*180/pi;
    while (phaseangle>180) || (abs(angle2vect2(av, a1)*180/pi) > ang_range(1))
        i = i+1;
        a1 = fcm(i, :);
        phaseangle = angle2vect(av, a1)*180/pi;
    end
    fcm(i, :) = [];
    dis(i) = [];
% removing all fcm along (0 degree) or against (180 degree) the a1 vector.
    ang = angle2vect(a1, fcm);
    ind = (abs(rem(ang*180/pi, 360)) < 10) | (abs(rem(ang*180/pi, 360)) >170);
    dis(ind) = [];
    fcm(ind, :) = [];

% finding b1

    % finding positive direction
    i = 1;
    b1 = fcm(i, :);
    [~, ind] = max(abs(b1));
    while (b1(ind) < 0) || (abs(angle2vect2(bv, b1)*180/pi) > ang_range(2))
        i = i+1;
        b1 = fcm(i, :);
        [~, ind] = max(abs(b1));
    end
    fcm(i, :) = [];
    dis(i) = [];
    
   
% finding c1.
% removing all fcm along (0 degree) or against (180 degree) the a1 vector.
    nv = -cross(a1, b1);
    ang = angle2vect(nv, fcm);
%    dis3 = dis(ind);
    % b1 should not be the second order of a1
    ind = (abs(rem(ang*180/pi, 360)) > 85) & (abs(rem(ang*180/pi, 360)) < 95);
    dis(ind) = [];
    fcm(ind,:) = [];
%     [~, ind] = min(dis);
%     c1 = fcm(ind, :);
    
    % finding positive direction
    i = 1;
    c1 = fcm(i, :);
    [~, ind] = max(abs(c1));
    while (c1(ind) < 0) || (abs(angle2vect2(cv, c1)*180/pi) > ang_range(3))
        i = i+1;
        c1 = fcm(i, :);
        [~, ind] = max(abs(c1));
    end
% Note that:
% [a1(:), b1(:), c1(:)] = latticevectormatrix' * [refa(:), refb(:), refc(:)]
% Therefore, 
% latticevectormat' = [a1(:),b1(:),c1(:)]*inv([refa(:), refb(:), refc(:)])
% latticevectormat = ([a1(:),b1(:),c1(:)]*inv([refa(:), refb(:),
% refc(:)]))'

Lmat = [a1(:),b1(:),c1(:)]/[refa(:), refb(:), refc(:)];
a1 = Lmat(:, 1)';
b1 = Lmat(:, 2)';
c1 = Lmat(:, 3)';
    set(handles.edview1, 'string', norm(a1));
    setappdata(gcbf, 'as', a1);
    c = 'r';
    linetag = 'linea';
    drawline(a1, linetag, c, handles);

    set(handles.edview2, 'string', norm(b1));
    setappdata(gcbf, 'bs', b1);
    c = 'g';
    linetag = 'lineb';
    drawline(b1, linetag, c, handles);

    set(handles.edview3, 'string', norm(c1));
    setappdata(gcbf, 'cs', c1);
    c = 'b';
    linetag = 'linec';
    drawline(c1, linetag, c, handles);
% construct a unit cell.
    pbprint_Callback([], [], handles);


% --------------------------------------------------------------------
function m_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to m_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ch_isReciprocal.
function ch_isReciprocal_Callback(hObject, eventdata, handles)
% hObject    handle to ch_isReciprocal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ch_isReciprocal
ev1 = findobj(gcbf, 'tag', 'edview1');
ev2 = findobj(gcbf, 'tag', 'edview2');
ev3 = findobj(gcbf, 'tag', 'edview3');
% alen = str2double(get(ev1, 'String'));
% blen = str2double(get(ev2, 'String'));
% clen = str2double(get(ev3, 'String'));

    if get(hObject,'Value')
        as = getappdata(handles.figure1, 'as');
        bs = getappdata(handles.figure1, 'bs');
        cs = getappdata(handles.figure1, 'cs');
    else
        as = getappdata(handles.figure1, 'a');
        bs = getappdata(handles.figure1, 'b');
        cs = getappdata(handles.figure1, 'c');
    end

    set(ev1, 'String', sprintf('%0.5f', norm(as)));
    set(ev2, 'String', sprintf('%0.5f', norm(bs)));
    set(ev3, 'String', sprintf('%0.5f', norm(cs)));
    

function v = get_vector_for_twopnts(varargin)
fh = getappdata(gcbf, 'dataFigure');
figure(fh);
ginput(1)
clickedpt1 = get(gca, 'CurrentPoint');
ginput(1)
clickedpt2 = get(gca, 'CurrentPoint');
%v1 = [clickedpt1(1, :),1]';
%v2 = VMtx*[clickedpt2(1, :),1]';
v1 = clickedpt1(1, :);
v2 = clickedpt2(1, :);
v = v2-v1;
vl = sqrt(sum(v.^2));
% make v normal to the vn
vn = get_view;
v = cross(cross(vn, v), vn);
v = v/norm(v)*vl;

% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v = get_vector_for_twopnts(handles);
set(handles.edview1, 'string', norm(v));
name = 'a';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
setappdata(gcbf, name, v);
%         setappdata(gcbf, 'aRmat', Rmat);
c = 'r';
linetag = 'linea';

drawline(v, linetag, c, handles);

% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v = get_vector_for_twopnts(handles);
set(handles.edview3, 'string', norm(v));
name = 'c';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
setappdata(gcbf, name, v);
%         setappdata(gcbf, 'aRmat', Rmat);
c = 'b';
linetag = 'linec';

drawline(v, linetag, c, handles);

% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v = get_vector_for_twopnts(handles);
set(handles.edview2, 'string', norm(v));
name = 'b';
if handles.ch_isReciprocal.Value
    name = [name, 's'];
end
setappdata(gcbf, name, v);
%         setappdata(gcbf, 'aRmat', Rmat);
c = 'g';
linetag = 'lineb';

drawline(v, linetag, c, handles);


% --------------------------------------------------------------------
function mn_integrateintensities_Callback(hObject, eventdata, handles)
% hObject    handle to mn_integrateintensities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = getappdata(gcbf, 'dataFigure');
hp = findobj(h, 'tag', 'grid_hkl');
hkls = get(hp, 'userdata');
d = sqrt(sum((hkls).^2, 2));
[~, ind] = min(d);

xd = hp.XData(:);
yd = hp.YData(:);
zd = hp.ZData(:);
pos = [xd(:), yd(:), zd(:)];
pos0 = pos(ind, :);
d = sqrt(sum((pos-pos0).^2, 2));
d(ind) = inf;
[~, ind] = min(d);
pos1 = pos(ind, :);

ud = get(h, 'userdata');
X = ud.Xax;
Y = ud.Yax;
Z = ud.Zax;
vox = ud.map;

qv = [X(:), Y(:), Z(:)];
d = sqrt(sum((qv-pos0).^2, 2));
[~, ind0] = min(d);

d = sqrt(sum((qv-pos1).^2, 2));
[~, ind1] = min(d);
[i0, j0, k0] = ind2sub(size(vox), ind0);
[i1, j1, k1] = ind2sub(size(vox), ind1);
ds = round(sqrt((i0-i1)^2+(j0-j1)^2+(k0-k1)^2)/2*0.75);

[inten, xp, yp, zp] = vox_hkl_integrate(vox, pos, X, Y, Z, ds);
% Fridel-pair averaging...
for i = 1:size(hkls, 1)
    G = hkls(i, 1:3);
    t = find(hkls(:,1) == -G(1) & hkls(:,2) == -G(2) & hkls(:,3) == -G(3));
    inten([i,t]) = sum(inten([i,t]))/2;
end
%end

ud.intensity = inten;
ud.qv = pos;
ud.xp = xp;
ud.yp = yp;
ud.zp = zp;
ud.hkls = hkls;
set(h, 'userdata', ud);
% % since qv = hkls*latticevectors', it may be possible to estimate
% % latticevectors from the fit.
% mat = ud.hkls'*[xp(:), yp(:), zp(:)]; % Gauss elimination...
% av = mat(:,1)';
% bv = mat(:,2)';
% cv = mat(:,3)';
% celcon4latvec(av, bv, cv)
options = optimset('fminsearch');
options = optimset(options, 'TolX',0.1E-8);
%        options = optimset(options, 'PlotFcns',@optimplotx);
options = optimset(options, 'MaxIter',5000);
options = optimset(options, 'MaxFunEvals', 5000);
cellinfo = getappdata(gcbf, 'cellinfo');
fprintf('Fitting optimum lattice vectors\n')
lv = fminsearch(@(x) fitf(x, [xp(:), yp(:), zp(:)], hkls),cellinfo.recilatticevectors, options);
fprintf('Fit result:')
cinf = celcon4rcplatvec(lv(1, :), lv(2, :), lv(3, :));
disp(cinf)
fprintf('Suggested lattice vectors ([av; bv; cv]):\n')
disp(cinf.latticevectors)
fprintf('Done.\n');

function err = fitf(x, qv, hkls)
% x is lattice vectors. [av; bv; cv];
qvc = hkls*x;
err = sum(sqrt(sum((qv-qvc).^2, 2)))/length(qvc);

% --------------------------------------------------------------------
function mn_saveintegratedintensities_Callback(hObject, eventdata, handles)
% hObject    handle to mn_saveintegratedintensities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = getappdata(gcbf, 'dataFigure');
ud = get(h, 'userdata');
if ~isfield(ud, 'hkls')
    fprintf('Run Analysis>Integrated Intensity before proceeding.\n');
end
data = [ud.hkls, ud.intensity/max(ud.intensity)*100, ones(size(ud.intensity))];
%sg = [];
filename = [];
cellinfo = getappdata(gcbf, 'cellinfo');

pmcobj = findobj(gcbf, 'tag', 'pm_centering_type');
str = get(pmcobj, 'String');
strv = get(pmcobj, 'Value');
[sgN, sg] = strtok(str{strv}, ':');
[sgN, syst] = strtok(sgN, ',');
sg = strtrim(sg(2:end));
syst = strtrim(syst(2:end));

if str2double(sgN)==0
    sg = [];
else
    switch syst
        case 'cubic'
            cellinfo.B = cellinfo.A;
            cellinfo.C = cellinfo.A;
            cellinfo.alpha = 90;
            cellinfo.beta = 90;
            cellinfo.gamma = 90;
        case {'hexagonal', 'trigonal'}
            cellinfo.B = cellinfo.A;
            cellinfo.alpha = 90;
            cellinfo.beta = 90;
            cellinfo.gamma = 120;
        case 'tetragonal'
            cellinfo.B = cellinfo.A;
            cellinfo.alpha = 90;
            cellinfo.beta = 90;
            cellinfo.gamma = 90;
        case 'orthorhombic'
            cellinfo.alpha = 90;
            cellinfo.beta = 90;
            cellinfo.gamma = 90;
        case 'monoclinic'
            cellinfo.alpha = 90;
            cellinfo.gamma = 90;
    end
end
saveasinflip([], data, sg, cellinfo, filename); % if data = [h,k,l,Iq,gN]


% --------------------------------------------------------------------
function mn_generatebkg_Callback(hObject, eventdata, handles)
% hObject    handle to mn_generatebkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = getappdata(gcbf, 'dataFigure');
hp = findobj(h, 'tag', 'grid_hkl');
hkls = get(hp, 'userdata');
d = sqrt(sum((hkls).^2, 2));
[~, ind] = min(d);

xd = hp.XData(:);
yd = hp.YData(:);
zd = hp.ZData(:);
pos = [xd(:), yd(:), zd(:)];
pos0 = pos(ind, :);
d = sqrt(sum((pos-pos0).^2, 2));
d(ind) = inf;
[~, ind] = min(d);
pos1 = pos(ind, :);

ud = get(h, 'userdata');
X = ud.Xax;
Y = ud.Yax;
Z = ud.Zax;
vox = ud.map;

qv = [X(:), Y(:), Z(:)];
d = sqrt(sum((qv-pos0).^2, 2));
[~, ind0] = min(d);

d = sqrt(sum((qv-pos1).^2, 2));
[~, ind1] = min(d);
[i0, j0, k0] = ind2sub(size(vox), ind0);
[i1, j1, k1] = ind2sub(size(vox), ind1);
ds = round(sqrt((i0-i1)^2+(j0-j1)^2+(k0-k1)^2)/2)+1;
vox = remove_intensity_from_vox(vox, pos, X, Y, Z, ds);
f = 1;
while f>0
    [vox, f] = fillNaN(vox);
end
ud.background = vox;
set(h, 'userdata', ud);

function calculate_q_pixel(varargin)
fh = getappdata(gcbf, 'dataFigure');
h = findobj(fh, 'tag', 'grid_hkl');
if isempty(h)
    return
end
if numel(h) == 2
    for i=1:numel(h)
        if ~isnan(h(i).MarkerFaceColor)
            h = h(i);
            break
        end
    end
end
xd = h.XData(:);
yd = h.YData(:);
zd = h.ZData(:);
qv = [xd, yd, zd];
%experimental setup
saxs = getappdata(gcbf, 'saxsinfo');
wl = saxs.waveln;
psize = saxs.psize;
sdd = saxs.SDD;
detangle = saxs.tiltangle;
phirng(1) = saxs.phirange(1);
phirng(2) = saxs.phirange(2);
center = saxs.center;
% calculating...
[pixN, phi] = qv2pixel(qv, wl, sdd, psize, detangle);
% plotting....
figh = getappdata(gcbf, 'expsumimgFigure');
ax = findobj(figh, 'type', 'axes');
xl = get(ax, 'XLim');
yl = get(ax, 'YLim');
delete(findobj(figh, 'tag', 'simulated_hkls'));
k = (phi>=phirng(1)) & (phi<=phirng(2)) & pixN(:,3) >= center(2);
k = k & pixN(:,2)+center(1)>=xl(1) & pixN(:,2)+center(1)<=xl(2);
k = k & pixN(:,3)+center(2)>=yl(1) & pixN(:,3)+center(2)<=yl(2);
h = plot(ax, pixN(k,2)+center(1), pixN(k,3)+center(2), 'ro');
k2 = ((phi<=(-180-phirng(1))) | (phi>=(180-phirng(2)))) & pixN(:,3) >= 0;
k2 = k2 & pixN(:,2)+center(1)>=xl(1) & pixN(:,2)+center(1)<=xl(2);
k2 = k2 & pixN(:,3)+center(2)>=yl(1) & pixN(:,3)+center(2)<=yl(2);
h2 = plot(ax, center(1)-pixN(k2,2), pixN(k2,3)+center(2), 'rs');
set(h, 'tag', 'simulated_hkls');
set(h2, 'tag', 'simulated_hkls');

% --------------------------------------------------------------------
function mn_indexing_Callback(hObject, eventdata, handles)
% hObject    handle to mn_indexing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mn_setcurrent_expsum_Callback(hObject, eventdata, handles)
% hObject    handle to mn_setcurrent_expsum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    img = evalin('base', 'img_sum');
    imagesc(log10(img+0.01), [5, 9]);axis image
    axis xy;
catch
    disp('img_sum is not available to plot. A current figure will be selected')
end
pbsetexpsumfigure_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function mn_overlay2D_Callback(hObject, eventdata, handles)
% hObject    handle to mn_overlay2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(hObject.Checked, 'on')
    set(hObject, 'Checked', 'Off');
else
    set(hObject, 'Checked', 'on');
    calculate_q_pixel
end


% --------------------------------------------------------------------
function mn_askexpsetup_Callback(hObject, eventdata, handles)
% hObject    handle to mn_askexpsetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% make a struct variable, mysetup, with the following fields
% mysetup
%   .waveln; for example 1 for 1A.
%   .psize; for example 0.075 for eiger
%   .SDD; for example 2000 for 2m
%   .tiltangle; for example = [0, 0, 0];
%   .phirange; for exasmple = [-80, 80];
try
    saxs = evalin('base', 'mysetup');
    setappdata(gcbf, 'saxsinfo', saxs);
catch
    disp('Make a struct variable, mysetup, with the following fields')
    disp(' mysetup')
    disp('   .waveln; for example 1 for 1A.')
    disp('   .psize; for example 0.075 for eiger')
    disp('   .SDD; for example 2000 for 2m')
    disp('   .tiltangle; for example = [0, 0, 0];')
    disp('   .phirange; for exasmple = [-80, 80];')
end


% --- Executes on button press in cb_show.
function cb_show_Callback(hObject, eventdata, handles)
% hObject    handle to cb_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_show
if get(hObject,'Value') 
    ct = getappdata(gcbf, 'cellinfo');
    v = eval(get(handles.edsetview, 'string'));
    v = v*ct.recilatticevectors;
    %setappdata(gcbf, 'as', as);
    c = 'r';
    linetag = 'line_hkl';
    h = drawline(v, linetag, c, handles);
    set(h, 'linewidth', 2)
    set(h, 'linestyle', ':')
else
    delete(findobj(0, 'tag', 'line_hkl'))
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider14_Callback(hObject, eventdata, handles)
% hObject    handle to slider14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider15_Callback(hObject, eventdata, handles)
% hObject    handle to slider15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider16_Callback(hObject, eventdata, handles)
% hObject    handle to slider16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider17_Callback(hObject, eventdata, handles)
% hObject    handle to slider17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ed_a_Callback(hObject, eventdata, handles)
% hObject    handle to ed_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_a as text
%        str2double(get(hObject,'String')) returns contents of ed_a as a double
update_to_cellinfo(handles)

% --- Executes during object creation, after setting all properties.
function ed_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_b_Callback(hObject, eventdata, handles)
% hObject    handle to ed_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_b as text
%        str2double(get(hObject,'String')) returns contents of ed_b as a double
update_to_cellinfo(handles)

% --- Executes during object creation, after setting all properties.
function ed_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_c_Callback(hObject, eventdata, handles)
% hObject    handle to ed_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_c as text
%        str2double(get(hObject,'String')) returns contents of ed_c as a double
update_to_cellinfo(handles)

% --- Executes during object creation, after setting all properties.
function ed_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to ed_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_alpha as text
%        str2double(get(hObject,'String')) returns contents of ed_alpha as a double
update_to_cellinfo(handles)

% --- Executes during object creation, after setting all properties.
function ed_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_beta_Callback(hObject, eventdata, handles)
% hObject    handle to ed_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_beta as text
%        str2double(get(hObject,'String')) returns contents of ed_beta as a double
update_to_cellinfo(handles)

% --- Executes during object creation, after setting all properties.
function ed_beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to ed_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_gamma as text
%        str2double(get(hObject,'String')) returns contents of ed_gamma as a double
update_to_cellinfo(handles)

% --- Executes during object creation, after setting all properties.
function ed_gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_psi_euler_Callback(hObject, eventdata, handles)
% hObject    handle to ed_psi_euler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_psi_euler as text
%        str2double(get(hObject,'String')) returns contents of ed_psi_euler as a double
update_to_cellinfo(handles)

% --- Executes during object creation, after setting all properties.
function ed_psi_euler_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_psi_euler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_hkl1.
function set_hkl1_Callback(hObject, eventdata, handles)
% hObject    handle to set_hkl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellinfo = getappdata(handles.figure1, 'cellinfo');
U = cellinfo.U;
hkl1 = eval(handles.ed_hkl1.String);
hkl = hkl1*cellinfo.recilatticevectors; % converting row reci vectors to cartesian. 
v = get_view; % get cartesian coordinates of the current view.
R = rotate_between_vectors(hkl/norm(hkl), v/norm(v));
newrv = cellinfo.recilatticevectors*R;
cellinfo = celcon4rcplatvec(newrv(1, :),newrv(2, :),newrv(3, :));
cellinfo.U = R'*U;
%cellinfo = orientcrystals(cellinfo, ori);
setappdata(handles.figure1, 'cellinfo', cellinfo);
update_from_cellinfo(handles)


function ed_hkl1_Callback(hObject, eventdata, handles)
% hObject    handle to ed_hkl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_hkl1 as text
%        str2double(get(hObject,'String')) returns contents of ed_hkl1 as a double


% --- Executes during object creation, after setting all properties.
function ed_hkl1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_hkl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_hkl2_Callback(hObject, eventdata, handles)
% hObject    handle to ed_hkl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_hkl2 as text
%        str2double(get(hObject,'String')) returns contents of ed_hkl2 as a double


% --- Executes during object creation, after setting all properties.
function ed_hkl2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_hkl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in set_hkl2.
function set_hkl2_Callback(hObject, eventdata, handles)
% hObject    handle to set_hkl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellinfo = getappdata(handles.figure1, 'cellinfo');
U = cellinfo.U;
hkl1 = eval(handles.ed_hkl1.String);
hkl2 = eval(handles.ed_hkl2.String);
hkl1 = hkl1*cellinfo.recilatticevectors; % converting row reci vectors to cartesian. 
hkl2 = hkl2*cellinfo.recilatticevectors; % converting row reci vectors to cartesian. 
hkl3 = cross(hkl1, hkl2);
v1 = hkl1;
v2 = get_view; % get cartesian coordinates of the current view.
v3 = cross(v1, v2);
R = [v1(:)/norm(v1), v2(:)/norm(v2), v3(:)/norm(v3)]*inv([hkl1(:)/norm(hkl1), hkl2(:)/norm(hkl2), hkl3(:)/norm(hkl3)]);
newrv = cellinfo.recilatticevectors*R';
newrv(1, :) = newrv(1, :)/norm(newrv(1, :))*norm(cellinfo.recilatticevectors(1, :)); 
newrv(2, :) = newrv(2, :)/norm(newrv(2, :))*norm(cellinfo.recilatticevectors(2, :)); 
newrv(3, :) = newrv(3, :)/norm(newrv(3, :))*norm(cellinfo.recilatticevectors(3, :)); 
cellinfo = celcon4rcplatvec(newrv(1, :),newrv(2, :),newrv(3, :));
%R = R*cellinfo.U;
cellinfo.U = (R*U);
setappdata(handles.figure1, 'cellinfo', cellinfo);
update_from_cellinfo(handles)


% --------------------------------------------------------------------
function load_cellinfo_frombase_Callback(hObject, eventdata, handles)
% hObject    handle to load_cellinfo_frombase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellinfo = evalin('base', 'cellinfo');
load_cellinfo(handles, cellinfo)
update_from_cellinfo(handles)


% --------------------------------------------------------------------
function exportCM2DataFigure_Callback(hObject, eventdata, handles)
% hObject    handle to exportCM2DataFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fh = getappdata(gcbf, 'dataFigure');
ud = get(fh, 'userdata');
udf = get(gcbf, 'userdata');
ud.CM = ud.CM+udf.cm/ud.pixelsize;
set(fh, 'userdata', ud);

% --------------------------------------------------------------------
function defineCM_Callback(hObject, eventdata, handles)
% hObject    handle to defineCM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fh = getappdata(gcbf, 'dataFigure');
ud = get(gcbf, 'userdata');
ax = findobj(fh, 'type', 'axes');
xl = ax.XLim;
yl = ax.YLim;
zl = ax.ZLim;
cm(1) = mean(xl);
cm(2) = mean(yl);
cm(3) = mean(zl);
ud.cm = cm;
set(gcbf, 'userdata', ud);


% --------------------------------------------------------------------
function definecenterofcoordinates_Callback(hObject, eventdata, handles)
% hObject    handle to definecenterofcoordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cb_show_cross.
function cb_show_cross_Callback(hObject, eventdata, handles)
% hObject    handle to cb_show_cross (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_show_cross
%h = findobj(0, 'tag', 'VolIqPlot');
%han = guihandles(h);
%fh = getappdata(gcbf, 'dataFigure');
if ~get(hObject, 'value')
    h = findobj(0, 'tag', 'VolIqPlot');
    delete(findobj(h, 'tag', 'hkl_slices'));
    return
end
switch lower(handles.bgview_cross.SelectedObject.String)
    case 'v2 x v3'
        v = [1, 0, 0];
    case 'v3 x v1'
        v = [0, 1, 0];
    case 'v1 x v2'
        v = [0, 0, 1];
end


handles.ed_tmp.Callback(hObject,v)
%draw_3dmap('drawhklslices',hObject,eventdata,handles)


function edsetview_cross_Callback(hObject, eventdata, handles)
% hObject    handle to edsetview_cross (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edsetview_cross as text
%        str2double(get(hObject,'String')) returns contents of edsetview_cross as a double


% --- Executes during object creation, after setting all properties.
function edsetview_cross_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edsetview_cross (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_tmp_Callback(hObject, eventdata, handles)
% hObject    handle to ed_tmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_tmp as text
%        str2double(get(hObject,'String')) returns contents of ed_tmp as a double


% --- Executes during object creation, after setting all properties.
function ed_tmp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_tmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function grid_q_max_Callback(hObject, eventdata, handles)
% hObject    handle to grid_q_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_q_max as text
%        str2double(get(hObject,'String')) returns contents of grid_q_max as a double


% --- Executes during object creation, after setting all properties.
function grid_q_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grid_q_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function mn_find_spacegroup_Callback(hObject, eventdata, handles)
% hObject    handle to mn_find_spacegroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
indexfound = find_allow_reflections(handles);
pmcobj = findobj(gcbf, 'tag', 'pm_centering_type');
str = get(pmcobj, 'String');
strv = get(pmcobj, 'Value');
% [N, restofN] = strtok(str{strv}, ',');
str_lattice = split(str, ', ');
N = str2double(str_lattice{strv,1});
sg_sel = split(str_lattice{strv,2}, ':');
latticetype = sg_sel{1};
for i=strv:size(str_lattice,1)
    sg_str = split(str_lattice{i,2}, ':');
    if ~strcmp(latticetype, sg_str{1})
        continue
    end
    sg = sg_gen(strtrim(sg_str{2}));
    isok = true;
    for k = size(indexfound, 1)
        if ccp4spg_is_sysabs(sg, indexfound(k, :))
            isok = false;
            break
        end
    end
    if isok
        fprintf('Possible spacegroup: %s\n', sg_str{2})
    end
end


% --------------------------------------------------------------------
function mn_plotobserved_reflections_Callback(hObject, eventdata, handles)
% hObject    handle to mn_plotobserved_reflections (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ff = evalin('base', 'QmapFamily');
if isempty(ff)
    warndlg('On the menu of draw_3dmap, run Tools>Split Isosurfaces.', 'Warning', 'modal');
    return
end

%fh = getappdata(gcbf, 'dataFigure');
%ax = findobj(fh, 'type', 'axes');
%ff = findobj(ax, 'type', 'patch');
fcm = zeros(numel(ff), 3);
for i=1:numel(ff)
    fcm(i, :) = ff{i}.cm;
end
Ndist = distNzero(fcm);
figure;
subplot(1,2,1)
plot3(fcm(:,1),fcm(:,2),fcm(:,3), 'ro');
title('QmapFamily cm')
axis image
subplot(1,2,2)
histogram(Ndist, round(numel(Ndist)/2))
title('Histogram of \Delta q')

% --------------------------------------------------------------------
function mn_swapaxis_Callback(hObject, eventdata, handles)
% hObject    handle to mn_swapaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mn_swapab_Callback(hObject, eventdata, handles)
% hObject    handle to mn_swapab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.ch_isReciprocal.Value
    reci_str = 's';
else
    reci_str = '';
end
ax = strcat('a', reci_str);
bx = strcat('b', reci_str);
cx = strcat('c', reci_str);
v1 = getappdata(gcbf, ax);
v2 = getappdata(gcbf, bx);
v3 = getappdata(gcbf, cx);
v = v2;
v2 = v1;
v1 = v;
setappdata(gcbf, ax, v1);
setappdata(gcbf, bx, v2);
setappdata(gcbf, cx, v3);

% --------------------------------------------------------------------
function mn_swapac_Callback(hObject, eventdata, handles)
% hObject    handle to mn_swapac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.ch_isReciprocal.Value
    reci_str = 's';
else
    reci_str = '';
end
ax = strcat('a', reci_str);
bx = strcat('b', reci_str);
cx = strcat('c', reci_str);
v1 = getappdata(gcbf, ax);
v2 = getappdata(gcbf, bx);
v3 = getappdata(gcbf, cx);
v = v3;
v3 = v1;
v1 = v;
setappdata(gcbf, ax, v1);
setappdata(gcbf, bx, v2);
setappdata(gcbf, cx, v3);


% --------------------------------------------------------------------
function mn_swapbc_Callback(hObject, eventdata, handles)
% hObject    handle to mn_swapbc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.ch_isReciprocal.Value
    reci_str = 's';
else
    reci_str = '';
end
ax = strcat('a', reci_str);
bx = strcat('b', reci_str);
cx = strcat('c', reci_str);
v1 = getappdata(gcbf, ax);
v2 = getappdata(gcbf, bx);
v3 = getappdata(gcbf, cx);
v = v2;
v2 = v3;
v3 = v;
setappdata(gcbf, ax, v1);
setappdata(gcbf, bx, v2);
setappdata(gcbf, cx, v3);
