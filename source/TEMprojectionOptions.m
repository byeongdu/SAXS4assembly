function varargout = TEMprojectionOptions(varargin)
% TEMPROJECTIONOPTIONS MATLAB code for TEMprojectionOptions.fig
%      TEMPROJECTIONOPTIONS, by itself, creates a new TEMPROJECTIONOPTIONS or raises the existing
%      singleton*.
%
%      H = TEMPROJECTIONOPTIONS returns the handle to a new TEMPROJECTIONOPTIONS or the handle to
%      the existing singleton*.
%
%      TEMPROJECTIONOPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEMPROJECTIONOPTIONS.M with the given input arguments.
%
%      TEMPROJECTIONOPTIONS('Property','Value',...) creates a new TEMPROJECTIONOPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TEMprojectionOptions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TEMprojectionOptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TEMprojectionOptions

% Last Modified by GUIDE v2.5 27-Mar-2025 16:23:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TEMprojectionOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @TEMprojectionOptions_OutputFcn, ...
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


% --- Executes just before TEMprojectionOptions is made visible.
function TEMprojectionOptions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TEMprojectionOptions (see VARARGIN)

% Choose default command line output for TEMprojectionOptions
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
if numel(varargin)>1
    if strcmp(get(varargin{1}, 'type'), 'uimenu')
        if strcmp(varargin{1}.Text, 'Film Geometry Tool')
            pbsetdatafigure_Callback(hObject, gcbf, handles)
        end
    end
end
% UIWAIT makes TEMprojectionOptions wait for user response (see UIRESUME)
% uiwait(handles.TEMprojection);


% --- Outputs from this function are returned to the command line.
function varargout = TEMprojectionOptions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function eddx_Callback(hObject, eventdata, handles)
% hObject    handle to eddx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eddx as text
%        str2double(get(hObject,'String')) returns contents of eddx as a double


% --- Executes during object creation, after setting all properties.
function eddx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eddx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edxf_Callback(hObject, eventdata, handles)
% hObject    handle to edxf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edxf as text
%        str2double(get(hObject,'String')) returns contents of edxf as a double


% --- Executes during object creation, after setting all properties.
function edxf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edxf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edyf_Callback(hObject, eventdata, handles)
% hObject    handle to edyf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edyf as text
%        str2double(get(hObject,'String')) returns contents of edyf as a double


% --- Executes during object creation, after setting all properties.
function edyf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edyf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edhklZ_Callback(hObject, eventdata, handles)
% hObject    handle to edhklZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edhklZ as text
%        str2double(get(hObject,'String')) returns contents of edhklZ as a double


% --- Executes during object creation, after setting all properties.
function edhklZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edhklZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edhklX_Callback(hObject, eventdata, handles)
% hObject    handle to edhklX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edhklX as text
%        str2double(get(hObject,'String')) returns contents of edhklX as a double


% --- Executes during object creation, after setting all properties.
function edhklX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edhklX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edmu_Callback(hObject, eventdata, handles)
% hObject    handle to edmu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edmu as text
%        str2double(get(hObject,'String')) returns contents of edmu as a double

if isappdata(gcbf, 'projectedImage')
    img = getappdata(gcbf, 'projectedImage');
    fobj = getappdata(gcbf, 'projectedImageFigureTag');
end
%fobj = findobj('tag', 'projectedImageFigureTag');
if ~isempty(fobj)
    mu = str2double(get(hObject,'String'));
    img = exp(-mu*img);
    imgh = findobj(fobj, 'type', 'image');
    set(imgh, 'CData', img);
end

% --- Executes during object creation, after setting all properties.
function edmu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edmu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edthick_Callback(hObject, eventdata, handles)
% hObject    handle to edthick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edthick as text
%        str2double(get(hObject,'String')) returns contents of edthick as a double


% --- Executes during object creation, after setting all properties.
function edthick_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edthick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eddz_Callback(hObject, eventdata, handles)
% hObject    handle to eddz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eddz as text
%        str2double(get(hObject,'String')) returns contents of eddz as a double


% --- Executes during object creation, after setting all properties.
function eddz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eddz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_projectTEM.
function pb_projectTEM_Callback(hObject, eventdata, handles)
% hObject    handle to pb_projectTEM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
projTEM(handles)

% --- Executes on button press in pb_make3D.
function pb_make3D_Callback(hObject, eventdata, handles)
% hObject    handle to pb_make3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
make3D(handles)


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edmu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edmu as text
%        str2double(get(hObject,'String')) returns contents of edmu as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edmu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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
    set(handles.pbsetdatafigure, 'visible', 'off');

function make3D(hd)
    fh = getappdata(gcbf, 'dataFigure');
    cellinfo = getappdata(fh, 'cellinfo');
    ud = get(fh, 'userdata');
    map3D = ud.map;
    Xax = ud.X;
    Yax = ud.Y;
    Zax = ud.Z;
    
    %fobj = findobj(0, 'tag', 'TEMprojection');
    %hd = guihandles(fobj);
    xf = str2double(hd.edxf.String);
    yf = str2double(hd.edyf.String);
    dx = str2double(hd.eddx.String);
    dz = str2double(hd.eddz.String);
    thick = str2double(hd.edthick.String);
    mu = str2double(hd.edmu.String);
    hklZ = eval(hd.edhklZ.String);
    hklX = eval(hd.edhklX.String);
    xx = 0:dx:xf;
    yy = 0:dx:yf;
    zz = 0:dz:thick;
    
    % dark
    IsoSurfLevel = getappdata(fh, 'isovalue');
    mp = map3D;
    if hd.rd_below.Value
        mp(map3D <= IsoSurfLevel) = 1;
        mp(map3D > IsoSurfLevel) = 0;
    else
        mp(map3D <= IsoSurfLevel) = 0;
        mp(map3D > IsoSurfLevel) = 1;
    end
    Z0 = hklZ;
    X0 = hklX;
    vtype = 'R = ';
    if hd.rb_real.Value
        hklZ = hklZ*cellinfo.mat';
        % % or Z = Z*cellinfo.latticevectors;
        hklX = hklX*cellinfo.mat';
        % % or Z = Z*cellinfo.latticevectors;
        vtype = 'R_{hkl} = ';
    end
    if hd.rb_rcp.Value
        hklZ = hklZ*cellinfo.recimat';
        % % or Z = Z*cellinfo.latticevectors;
        hklX = hklX*cellinfo.recimat';
        % % or Z = Z*cellinfo.latticevectors;
        vtype = 'G_{hkl} = ';
    end
    
    isovalue = getappdata(fh, 'isovalue');
    [Vox, fm] = voxreconstruct(xx, yy, zz, hklZ, hklX, {Xax, Yax, Zax, map3D}, ...
        cellinfo, isovalue, 'absolute');
    xlabel(sprintf('%s[%i,%i,%i] direction',vtype, X0), 'fontsize', 15, 'color', 'g')
    zlabel(sprintf('%s[%i,%i,%i] direction',vtype, Z0), 'fontsize', 15, 'color', 'r')
    U = rotate_plane_onto_XY(hklZ, hklX);
    dz = dspacing(Z0, cellinfo, 'triclinic');
    %dz = dspacing('triclinic', Z0, cellinfo);
    dx = dspacing(X0, cellinfo, 'triclinic');
    hklZ = hklZ*inv(U);
    hklz = hklZ/norm(hklZ)*dz;
    hklx = hklX*inv(U)/norm(hklX)*dx;
    O = [0, 0, 0];
    k=arrow3(O, hklz, 'r1', 0.2);
    set(k(1), 'linestyle', '-', 'color', 'r');
    set(k, 'tag', 'hklvector');
    k=arrow3(O, hklx, 'g1', 0.2);
    set(k(1), 'linestyle', '-', 'color', 'g');
    set(k, 'tag', 'hklvector');
    dt = [];
    dt.isovalue = isovalue;
    dt.U = U;
    dt.zaxis = zz;
    dt.Vox = Vox;
    dt.cellinfo = cellinfo;
    set(gcf, 'userdata', dt);
    f2 = uimenu(fm, 'Label','Process');
    uimenu('Parent',f2,'Label','Convert Patches to Objects','Callback',@Patch2Object);
    uimenu('Parent',f2,'Label','Integrate XY planes','Callback',@integrateXYplane);
    function integrateXYplane(varargin)
        m = findobj(gcf, 'type', 'patch');
        dt = get(gcbf, 'userdata');
        Vox = dt.Vox;
        cgy = sum(Vox, 1);
        cgy = sum(cgy, 2);
        cgy = squeeze(cgy);
        figure;
        plot(dt.zaxis, cgy/max(cgy), 'r');
        
    function Patch2Object(varargin)
        m = findobj(gcf, 'type', 'patch');
        dt = get(gcbf, 'userdata');
        [face, verts] = convertpatch2object(m);
        [~,I] = sort(cellfun(@length,face),'descend');
        face = face(I);
        Nel = cellfun(@length,face);
        Ncolor = sum(Nel>mean(Nel));
        h = figure;
        cols = jet(Ncolor);
        for i=1:numel(face)
            if i<=Ncolor
                c = cols(i, :);
            else
                c = cols(Ncolor, :);
            end
            patch('faces',face{i},...
                        'vertices', verts,...
                        'FaceColor',c, 'edgecolor', 'none');
        end
        axis image;
        camlight;
        lighting gouraud;
        camva(8);
        hold on;
        setappdata(h, 'cellinfo', dt.cellinfo);
        setappdata(h, 'RotMat', dt.U);
        f2 = uimenu(h, 'Label','View');
        uimenu('Parent',f2,'Label','Draw Unitcell','Callback',@drawcell);
        uimenu('Parent',f2,'Label','HKL Drawing Tool','Callback',@hkldrawingtool);
        uimenu('Parent',f2,'Label','Color Objects','Callback',{@setcolor, h});
    function setcolor(varargin)
        colorobjects(varargin{3});
    
    function drawcell(varargin)
        ci = getappdata(gcbf, 'cellinfo');
        mat = ci.mat;
        try
            U = getappdata(gcbf, 'RotMat');
        catch
            U = eye(3);
        end
        lv = ci.latticevectors*inv(U);
        mat = lv';
        %cellinfo = celcon4latvec(lv(1, :), lv(2, :), lv(3, :));
        %mat = cellinfo.mat;
        %setappdata(gcbf, 'cellinfo', cellinfo);
%         if ~isempty(U)
%             mat = mat*U;
%         end
        switch get(varargin{1}, 'Checked')
            case 'off'
                a1 = mat(1, :);
                a2 = mat(2, :);
                a3 = mat(3, :);
                t = drawunitcell(a1, a2, a3);
                set(varargin{1}, 'Checked', 'on')
                set(t, 'Tag', 'mycell')
            case 'on'
                projh = findobj(gcbf, 'tag', 'mycell');
                delete(projh);
                set(varargin{1}, 'Checked', 'off')
        end        
    
function projTEM(hd)
    fh = getappdata(gcbf, 'dataFigure');
    cellinfo = getappdata(fh, 'cellinfo');
    ud = get(fh, 'userdata');
    map3D = ud.map;
    Xax = ud.X;
    Yax = ud.Y;
    Zax = ud.Z;
    
    %fobj = findobj(0, 'tag', 'TEMprojection');
    %hd = guihandles(fobj);
    xf = str2double(hd.edxf.String);
    yf = str2double(hd.edyf.String);
    dx = str2double(hd.eddx.String);
    dz = str2double(hd.eddz.String);
    thick = str2double(hd.edthick.String);
    mu = str2double(hd.edmu.String);
    hklZ = eval(hd.edhklZ.String);
    hklX = eval(hd.edhklX.String);
    xx = 0:dx:xf;
    yy = 0:dx:yf;
    % dark
    IsoSurfLevel = getappdata(fh, 'isovalue');
    mp = map3D;
    if hd.rd_below.Value
        mp(map3D <= IsoSurfLevel) = 1;
        mp(map3D > IsoSurfLevel) = 0;
    else
        mp(map3D <= IsoSurfLevel) = 0;
        mp(map3D > IsoSurfLevel) = 1;
    end
    
    Z0 = hklZ;
    X0 = hklX;

    vtype = 'R = ';
    if hd.rb_real.Value
        hklZ = hklZ*cellinfo.mat';
        % % or Z = Z*cellinfo.latticevectors;
        hklX = hklX*cellinfo.mat';
        % % or Z = Z*cellinfo.latticevectors;
        vtype = 'R_{hkl} = ';
    end
    if hd.rb_rcp.Value
        hklZ = hklZ*cellinfo.recimat';
        % % or Z = Z*cellinfo.latticevectors;
        hklX = hklX*cellinfo.recimat';
        % % or Z = Z*cellinfo.latticevectors;
        vtype = 'G_{hkl} = ';
    end
        
    img = projVoxel2Image({Xax, Yax, Zax, mp}, cellinfo, xx, yy, thick, hklZ, hklX, dz);
%        img = projVoxel2Image(mp, cellinfo, xx, yy, thick, hklZ, hklX, dz);
    TEMimg = exp(-mu*img);
    ftem = figure;
    set(ftem, 'tag', 'projectedTEMimage');
    imagesc(xx, yy, TEMimg);
    axis image; colorbar;
    axis xy;
    xlabel(sprintf('%s[%i,%i,%i] direction',vtype, X0), 'fontsize', 15)
    title(sprintf('Looking through %s[%i,%i,%i] direction', vtype, Z0), 'fontsize', 15)
    setappdata(gcbf, 'projectedImage', img);
    setappdata(gcbf, 'projectedImageFigureTag', ftem);
 


% --- Executes on button press in pb_make1repeating.
function pb_make1repeating_Callback(hObject, eventdata, handles)
% hObject    handle to pb_make1repeating (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fh = getappdata(gcbf, 'dataFigure');
    cellinfo = evalin('base', 'cellinfo');
    ud = get(fh, 'userdata');
    
        
    %fobj = findobj(0, 'tag', 'TEMprojection');
    %hd = guihandles(fobj);
    hklZ = eval(handles.edhklZ.String);
    hklX = eval(handles.edhklX.String);
    [Vq, hklx, hkly, hklz, cx, cy, cz] = transformUC(ud.X, ud.Y, ud.Z, ud.map,hklZ, cellinfo, hklX);

    vtype = 'G_{hkl} = ';
    
    h = figure;
    x = linspace(0, 1, size(Vq, 1)+1)*cx;
    y = linspace(0, 1, size(Vq, 2)+1)*cy;
    z = linspace(0, 1, size(Vq, 3)+1)*cz;
    x(1) = [];
    y(1) = [];
    z(1) = [];
    [Xx, Yy, Zz] = ndgrid(x, y, z);
    isovalue = mean(Vq(:));
    [faces, verts, cols] = isosurface(Xx, Yy, Zz, Vq, isovalue, Zz);
    [F,V,C] = isocaps(Xx, Yy, Zz, Vq, isovalue);
    llmap = min(Vq(:));
    ulmap = max(Vq(:));
    cols = (isovalue-llmap)/(ulmap-llmap)*255+1;
    cmp = jet(256);
    cols = cmp(round(cols), :);

    % channels..
    patch('faces',faces,...
                    'vertices', verts,...
                    'FaceColor',cols, 'edgecolor', 'none');
    % patch('faces',faces,...
    %                 'vertices', verts,...
    %                 'facevertexCData',cols,...
    %                 'FaceColor','interp', 'edgecolor', 'none');
    hold on
    % Caps...

    cols = (C-llmap)/(ulmap-llmap)*255+1;
    cols = cmp(round(cols), :);
    patch('faces',F,...
                    'vertices', V,...
                    'facevertexCData',cols,...
                    'FaceColor','interp', 'edgecolor', 'none');
    %ylabel(sprintf('Ghkl=[%i,%i,%i] direction (nm)',Ydir), 'fontsize', 15)
    %zlabel(sprintf('Ghkl=[%i,%i,%i] direction (nm)',Zdir), 'fontsize', 15)
    
    axis image;
    camlight;
    lighting gouraud;
    camva(8);

    xlabel(sprintf('<%i,%i,%i>',hklx), 'fontsize', 15, 'color', 'g')
    zlabel(sprintf('<%i,%i,%i>',hklz), 'fontsize', 15, 'color', 'r')
    ylabel(sprintf('<%i,%i,%i>',hkly), 'fontsize', 15, 'color', 'b')
    

    O = [0, 0, 0];
    k=arrow3(O, hklz, 'r1', 0.2);
    set(k(1), 'linestyle', '-', 'color', 'r');
    set(k, 'tag', 'hklvector');
    k=arrow3(O, hklx, 'g1', 0.2);
    set(k(1), 'linestyle', '-', 'color', 'g');
    set(k, 'tag', 'hklvector');
    dt = [];
    dt.isovalue = isovalue;
    dt.zaxis = z;
    dt.Vox = Vq;
    dt.cellinfo = cellinfo;
    set(gcf, 'userdata', dt);
    f2 = uimenu(h, 'Label','Process');
    uimenu('Parent',f2,'Label','Convert Patches to Objects','Callback',@Patch2Object);
    uimenu('Parent',f2,'Label','Integrate XY planes','Callback',@integrateXYplane);