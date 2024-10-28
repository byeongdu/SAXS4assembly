function setgihandle(varargin)
% This will save 'saxs' on 'SAXSImageViewer' as well as 'GISAXSLee'
if numel(varargin) >= 1
    saxs = varargin{1};
    handle = [];
end

if numel(varargin) == 2
    handle = varargin{2};
end

try
    hSAXSlee = evalin('base', 'SAXSimageviewerhandle');
catch
    disp('SAXSimageviewer does not exist....');
    hSAXSlee = [];
end

%if isempty(handle)
%    hSAXSlee=findobj('Tag','SAXSImageViewer');
%else
%    hSAXSlee=findobj(handle, 'Tag','SAXSImageViewer');
%end

if ~isempty(hSAXSlee)
    set(hSAXSlee, 'Userdata', saxs);
else
    set(gcbf, 'Userdata', saxs);    
end

hSAXSlee=findobj('Tag','gisaxsleenew');

if ~isempty(hSAXSlee)
    set(hSAXSlee, 'Userdata', saxs);
end

%if isfield(saxs, 'imgfigurehandle')
%    if ~isempty(saxs.imgfigurehandle)
%        try
%            set(saxs.imgfigurehandle, 'Userdata', saxs);
%        catch
%            strm = sprintf('there is no figure handle %i in saxs\n', saxs.imgfigurehandle);
%            strm = sprintf('%sError: in setgihandle.m', strm);
%            disp(strm)
%        end
%    end
%end
