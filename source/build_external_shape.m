function varargout = build_external_shape(varargin)
% BUILD_EXTERNAL_SHAPE MATLAB code for build_external_shape.fig
%      BUILD_EXTERNAL_SHAPE, by itself, creates a new BUILD_EXTERNAL_SHAPE or raises the existing
%      singleton*.
%
%      H = BUILD_EXTERNAL_SHAPE returns the handle to a new BUILD_EXTERNAL_SHAPE or the handle to
%      the existing singleton*.
%
%      BUILD_EXTERNAL_SHAPE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BUILD_EXTERNAL_SHAPE.M with the given input arguments.
%
%      BUILD_EXTERNAL_SHAPE('Property','Value',...) creates a new BUILD_EXTERNAL_SHAPE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before build_external_shape_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to build_external_shape_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help build_external_shape

% Last Modified by GUIDE v2.5 05-Mar-2022 18:44:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @build_external_shape_OpeningFcn, ...
                   'gui_OutputFcn',  @build_external_shape_OutputFcn, ...
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


% --- Executes just before build_external_shape is made visible.
function build_external_shape_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to build_external_shape (see VARARGIN)

% Choose default command line output for build_external_shape
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes build_external_shape wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = build_external_shape_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ed_normalvector_Callback(hObject, eventdata, handles)
% hObject    handle to ed_normalvector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_normalvector as text
%        str2double(get(hObject,'String')) returns contents of ed_normalvector as a double
v = eval(sprintf('[%s]', hObject.String));
ud = get(gcf, 'userdata');
if isfield(ud, 'cellinfo')
    cellinfo = ud.cellinfo;
else
    cellinfo = evalin('base', 'cellinfo');
    fprintf('Cellinfo is loaded from workspace.\n');
end
%[D, HKLs] = indexref_hkl(cellinfo, sgroup('P 1'), v);
% for i=1:numel(HKLs)
%     hkls{i} = mat2str(HKLs(i).HKLs');
% end
hkl = get_allpermute_triplets(v(1), v(2), v(3));
for i=1:size(hkl,1)
    hkls{i} = mat2str(hkl(i, :));
end
set(handles.pm_hkls, 'string', hkls);
set(handles.pm_hkls, 'value', 1);

% --- Executes during object creation, after setting all properties.
function ed_normalvector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_normalvector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_add.
function pb_add_Callback(hObject, eventdata, handles)
% hObject    handle to pb_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v = eval(sprintf('[%s]', handles.ed_normalvector.String));
draw(v)

function draw(v, ishkl, P)
    if nargin==1
        ishkl = true;
    end
    if nargin<3
        P = [];
    end
    h = findobj(gcf, 'tag', 'currentplane');
    if ~isempty(h)
        delete(h)
    end
    cellinfo = evalin('base', 'cellinfo');
    %v = v*cellinfo.mat;
    if ishkl
        v = v*cellinfo.recilatticevectors; % real space vector
    end
    
    %v = [v(2), v(1), v(3)];
    box = [xlim, ylim, zlim];
    [O, dmax] = get_traveling_plane_in_box(v, box);
    if ~isempty(P)
        planeO = createPlane(P, v);
        d = distancePointPlane(O, planeO);
        f_d = abs(d)/dmax;
    else
        f_d = 0.5;
    end
    [pts, P] = traveling_plane_in_box(v, O, f_d, dmax, box);
    set(findobj(gcbf, 'tag', 'planeslider'), 'value', f_d);
    %h = drawPlane3d(pts);
    h = drawPolyhedron(pts, 1:size(pts, 1));
    set(h, 'facecolor', 'y', 'facealpha', 0.1, 'edgecolor', 'k');
    set(h, 'tag', 'currentplane')
    ud.O = O;
    ud.dmax = dmax;
    ud.box = box;
    ud.v = v;
    ud.P0 = P;
    set(h, 'userdata', ud);


% --- Executes on slider movement.
function planeslider_Callback(hObject, eventdata, handles)
% hObject    handle to planeslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    %v = eval(sprintf('[%s]', handles.ed_normalvector.String));
    h = findobj(gcf, 'tag', 'currentplane');
    if isempty(h)
        return
    end
    ud = get(h, 'userdata');
    O = ud.O;
    v = ud.v;
    dmax = ud.dmax;
    box = ud.box;
    pts = traveling_plane_in_box(v, O, hObject.Value, dmax, box);
    delete(h)
    %h = drawPlane3d(pts, 'facecolor', 'y');
    h = drawPolyhedron(pts, 1:size(pts, 1));
    set(h, 'facecolor', 'y', 'facealpha', 0.1, 'edgecolor', 'k');
    set(h, 'tag', 'currentplane')
    ud.O = O;
    ud.dmax = dmax;
    ud.box = box;
    set(h, 'userdata', ud);

% --- Executes during object creation, after setting all properties.
function planeslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to planeslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pb_remove.
function pb_remove_Callback(hObject, eventdata, handles)
% hObject    handle to pb_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
planes = get(handles.pm_planes, 'string');
val = get(handles.pm_planes, 'value');
%name = planes{val};
if iscell(planes)
    name = planes{val};
    planes(val) = [];
else
    name = planes;
    planes = '';
end
h = findobj(gcf, 'tag', name);
delete(h)
set(handles.pm_planes, 'string', planes);

% --- Executes on selection change in pm_planes.
function pm_planes_Callback(hObject, eventdata, handles)
% hObject    handle to pm_planes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_planes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_planes
planes = get(handles.pm_planes, 'string');
val = get(handles.pm_planes, 'value');
if iscell(planes)
    name = planes{val};
    planes(val) = [];
else
    name = planes;
    planes = '';
end
h = findobj(gcf, 'tag', name);
set(h, 'tag', 'currentplane');
set(handles.pm_planes, 'string', planes);

% --- Executes during object creation, after setting all properties.
function pm_planes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_planes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_keep.
function pb_keep_Callback(hObject, eventdata, handles)
% hObject    handle to pb_keep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
planes = get(handles.pm_planes, 'string');
h = findobj(gcf, 'tag', 'currentplane');
for i=1:1000
    name = sprintf('plane%i', i);
    if ~findmystr(planes, name)
        if iscell(planes)
            planes = [planes; name];
        else
            planes = {planes; name};
        end
        set(h, 'tag', name);
        set(h, 'facecolor', 'g', 'facealpha', 0.05)
        set(handles.pm_planes, 'string', planes);
        break
    end
end

function val = findmystr(celln, mystr)
    if iscell(celln)
        nc = numel(celln);
    else
        celln = {celln};
        nc = 1;
    end
    for i=1:nc
        if strcmp(celln{i}, mystr)
            val = true;
            return
        end
    end
    val = false;
    
% --- Executes on button press in pb_buildpolyhedron.
function pb_buildpolyhedron_Callback(hObject, eventdata, handles)
% hObject    handle to pb_buildpolyhedron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
planes = get(handles.pm_planes, 'string');
A = [];
b = [];
eqn = [];
for i=1:numel(planes)
    if isempty(planes{i})
        continue
    end
    obj = findobj(gcf, 'tag', planes{i});
    ud = get(obj, 'userdata');
    
    if isempty(ud)
        continue
    end
    
    [v, b0] = get_A_b(obj);
    
    A = [A;v];
    % [0, 0, 0] should be included..
    % since, a*x+b*y+c*z <= a*x0+b*y0+c*z0, b0 should be positive.
    b = [b; abs(b0)];
    delete(obj)
end
% Vortices of the polyhedron is confined by Ax <= b.
V=MY_con2vert(A,b);
%V=MY_con2vert(eqn(:, 1:3),eqn(:, end));
%V = [V(:,2), V(:,1), V(:,3)];
ind = convhull(V(:,1), V(:,2), V(:,3));
vert = V(ind, :);
vert = unique_m(vert);
%[V,nr] = MY_con2vert(A,b);
F = findfaces(vert);
alpha.vertices = vert;
alpha.faces = F;
alpha.color = 'r';
alpha.facealpha = 0.2;
assignin('base', 'extrn', alpha);
t = drawpolyhedron(alpha);
set(t, 'tag', 'ExternalShape');
set(handles.pm_planes, 'string', {});

function [v, b0] = get_A_b(varargin)
    obj = varargin{1};
    ud = get(obj, 'userdata');
    if isempty(ud)
        v = []; b0 = [];
        return
    end
    v = ud.v;
    center = [mean(ud.box(1:2)),mean(ud.box(3:4)),mean(ud.box(5:6))];
    % ud.v and ud.box
    if isfield(ud, 'P0')
        P0 = ud.P0;
    else
        P0 = mean(obj.Vertices);
    end
    
    % if the normal vector is pointing inward, the sign needs to be flipped. 
    if (P0-center)*v' < 0
        v = -v;
    end

    b0 = v*P0';

% --- Executes on button press in tg_drawunitcell.
function tg_drawunitcell_Callback(hObject, eventdata, handles)
% hObject    handle to tg_drawunitcell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value
    cellinfo = evalin('base', 'cellinfo');
    as = cellinfo.latticevectors(1, :);
    bs = cellinfo.latticevectors(2, :);
    cs = cellinfo.latticevectors(3, :);
    xl = xlim;
    yl = ylim;
    zl = zlim;
    O = [mean(xl), mean(yl), mean(zl)];
    scale = min([diff(xl), diff(yl), diff(zl)])/2;
    unc = max([norm(as), norm(bs), norm(cs)]);
    as = as/unc*scale;
    bs = bs/unc*scale;
    cs = cs/unc*scale;
    t = drawunitcell(as, bs, cs);
    set(t, 'Tag', 'mycell')
else
    delete(findobj(gcf, 'tag', 'mycell'))
end


% --- Executes on selection change in pm_hkls.
function pm_hkls_Callback(hObject, eventdata, handles)
% hObject    handle to pm_hkls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_hkls contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_hkls

contents = cellstr(get(hObject,'String'));
v = contents{get(hObject,'Value')};
draw(eval(v));

% --- Executes during object creation, after setting all properties.
function pm_hkls_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_hkls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_fromdatatip.
function pb_fromdatatip_Callback(hObject, eventdata, handles)
% hObject    handle to pb_fromdatatip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dt = findall(gcf, 'Type', 'DataTip');
pnts=[];
for i=1:numel(dt)
    pnts(i, :) = [dt(i).X, dt(i).Y, dt(i).Z];
end
%plane = createPlane(pnts);
plane = [pnts(1, :), pnts(2, :)-pnts(1, :), pnts(3, :)-pnts(1, :)];
v = planeNormal(plane);
draw(v, false);
delete(dt);


% --- Executes on button press in pb_addview.
function pb_addview_Callback(hObject, eventdata, handles)
% hObject    handle to pb_addview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
VMtx = view(gca);
ginput(1)
clickedpt1 = get(gca, 'CurrentPoint');
ginput(1)
clickedpt2 = get(gca, 'CurrentPoint');
%v1 = [clickedpt1(1, :),1]';
%v2 = VMtx*[clickedpt2(1, :),1]';
v1 = clickedpt1(1, :);
v2 = clickedpt2(1, :);
v = v2-v1;
[az, el] = view(gca);
%Rmat = viewmtx(az, el);
vz = sin(el*pi/180);
vx = cos(el*pi/180)*sin(az*pi/180);
vy = -cos(el*pi/180)*cos(az*pi/180);
as = [vx, vy, vz];%view(viewvector)

N = cross(as, v(1:3));
draw(N, false, v1);


% --------------------------------------------------------------------
function mn_File_Callback(hObject, eventdata, handles)
% hObject    handle to mn_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mn_loadcellinfo_Callback(hObject, eventdata, handles)
% hObject    handle to mn_loadcellinfo (see GCBO)
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
end
ud = get(gcf, 'userdata');
ud.cellinfo = cellinfo;
set(gcf, 'userdata', ud);
% --------------------------------------------------------------------
function mn_savepolyhedron_Callback(hObject, eventdata, handles)
% hObject    handle to mn_savepolyhedron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, pname] = uiputfile('*.mat', 'Save cellinfo As');
if isequal(fname,0) || isequal(pname,0)
   disp('User pressed cancel')
   return
end
fn = fullfile(pname, fname);
externalshape = evalin('base', 'extrn');
save(fn, 'externalshape')
fprintf('Saved as %s.\n', fn)


% --------------------------------------------------------------------
function mn_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to mn_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mn_remove_voxel_outsideplane_Callback(hObject, eventdata, handles)
% hObject    handle to mn_remove_voxel_outsideplane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~strcmp(get(gcf, 'Tag'), 'VolDataPlot')
    error('Current figure is not voxelmap')
end
ud = get(gcf, 'userdata');
obj = findobj(gcf, 'tag', 'currentplane');
[v, b0] = get_A_b(obj);
pos = [ud.Y(:)-ud.CM(2), ud.X(:)-ud.CM(1), ud.Z(:)-ud.CM(3)]*ud.pixelsize;
t = pos*v'>b0;
ud.map(t) = 0;
set(gcf, 'userdata', ud);
% --------------------------------------------------------------------
function mn_remove_voxel_belowvalue_Callback(hObject, eventdata, handles)
% hObject    handle to mn_remove_voxel_belowvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~strcmp(get(gcf, 'Tag'), 'VolDataPlot')
    error('Current figure is not voxelmap')
end
guh = guihandles(gcf);
ud = get(gcf, 'userdata');
t = ud.map < guh.isosurfslider.Value;
ud.map(t) = 0;
set(gcf, 'userdata', ud);