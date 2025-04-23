    function colorobjects(varargin)
        if numel(varargin) == 0
            fh = gcf;
        else
            fh = varargin{1};
        end
        m = findobj(fh, 'type', 'patch');
        npatches = linspace(1, numel(m), numel(m));
        
        %d = dialog('Position',[300 300 250 250],'Name','Select Color');
        d = figure('Position',[300 300 250 250],'Name','Select Color');
        set(d, 'menubar', 'none');
        set(d, 'toolbar', 'none');
    
        set(d, 'userdata', fh);
        uicontrol('Parent',d,...
               'Style','text',...
               'Position',[20 200 210 25],...
               'String','Choose Object(s), otherwise from GUI:');

        mn_sel = uicontrol('Parent',d,...
               'Style','popup',...
               'tag', 'pm_object',...
               'Position',[75 180 100 25],...
               'String',npatches,...
               'value', 1,...
               'Callback',@popup_chooseobject_callback);
           
        uicontrol('Parent',d,...
               'Style','text',...
               'tag', 'txt_colortext',...
               'Position',[20 140 210 25],...
               'String','Color for Selected Object(s):');
           
        uicontrol('Parent',d,...
               'Style','pushbutton',...
               'Position',[75 120 100 25],...
               'String','Select Color',...
               'Callback',@popup_choosecolor_callback);
           
        uicontrol('Parent',d,...
               'Style','text',...
               'Position',[20 80 210 25],...
               'String','Facealpha for Selected Object(s) :');
        uicontrol('Parent',d,...
               'Style','edit',...
               'tag','ed_facealpha',...
               'Position',[75 60 100 25],...
               'String','1',...
               'Callback',@popup_choosefacealpha_callback);
          
           
        uicontrol('Parent',d,...
               'Position',[89 20 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');
        
       %sel = findobj(fh, 'tag', 'objects', 'selected', 'on');
       v = [m.Selected];
       sel = find(v==1);
       if ~isempty(sel)
           mn_sel.Value = sel;
       end

        % Wait for d to close before running to completion
%        uiwait(d);
           function popup_chooseobject_callback(varargin)
               if ishandle(varargin{1})
                   popup = varargin{1};
                   idx = popup.Value;
                   fhd = get(gcbf, 'userdata');
                   oj = findobj(fhd, 'type', 'patch');
                   cl = get(oj(idx), 'facecolor');
                   fa = get(oj(idx), 'facealpha');
                   set(oj(idx), 'Selected', 'On');
                   set(findobj(gcbf, 'tag', 'txt_colortext'), ...
                       'foregroundcolor', cl);
                   set(findobj(gcbf, 'tag', 'ed_facealpha'), ...
                       'string', num2str(fa));
               end
           end
           function popup_choosecolor_callback(varargin)
               c = uisetcolor;
               %popup = findobj(gcbf, 'tag', 'pm_object');
               %idx = popup.Value;

               fhd =get(gcbf, 'userdata');
               oj = findobj(fhd, 'type', 'patch', 'selected', 'on');
               %set(oj(idx), 'facecolor', c);
               set(oj, 'facecolor', c);
               set(findobj(gcbf, 'tag', 'txt_colortext'), ...
                   'foregroundcolor', c);
               set(oj, 'Selected', 'off')
           end
           function popup_choosefacealpha_callback(varargin)
               %popup = findobj(gcbf, 'tag', 'pm_object');
               %idx = popup.Value;
               %fhd =get(gcbf, 'userdata');
               fhd =get(gcbf, 'userdata');
               oj = findobj(fhd, 'type', 'patch', 'selected', 'on');
               set(oj, 'facealpha', str2double(varargin{1}.String));
               set(oj, 'Selected', 'off')
               %set(oj(idx), 'facealpha', str2double(varargin{1}.String));
%               set(oj(idx), 'Selected', 'Off');
           end
    end
