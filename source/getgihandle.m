function saxs = getgihandle(varargin)
if isempty(varargin)
    handle = [];
else
    handle = varargin{1};
end

try 
    imgh = evalin('base', 'SAXSimageviewerhandle');
catch
    disp('SAXSimageviewer does not exist....');
    imgh = [];
end

%if isempty(handle)
%    imgh = findobj('Tag', 'SAXSImageViewer');
%else
%    imgh = findobj(handle, 'Tag', 'SAXSImageViewer');
%end

% it will first try to read 'saxs' from 'SAXSImageViewer'
if ~isempty(imgh)
    saxs = get(imgh, 'Userdata');
    % if SAXSImageViewer has 'saxs', it will stop reading saxs in GISAXSLee.
    if ~isempty(saxs)
        return
    end
end

% When SAXSImageViewer is empty, it will try GISAXSLee
hSAXSlee = findobj('Tag','gisaxsleenew');
saxs=get(hSAXSlee, 'Userdata');

end