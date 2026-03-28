function add_imagemenu(varargin)
% add_imagemenu  Add an "Image" menu to a figure containing an imagesc plot.
%
% Usage: add_imagemenu(fig)
%
% Before calling, store the raw linear-scale data and current mode:
%   setappdata(fig, 'imgraw',   linear_data);
%   setappdata(fig, 'imgislog', true/false);   % true = currently log scale

menuTag = 'figMenuImage';
if ~isempty(varargin)
    fig = varargin{1};
    if isempty(findobj(fig, 'tag', menuTag))
        f = uimenu(fig, 'Label', 'Image');
        set(f, 'tag', menuTag)
    else
        return
    end
else
    fig = gcf;
    if isempty(findobj(fig, 'tag', menuTag))
        f = uimenu('Label', 'Image');
        set(f, 'tag', menuTag)
    else
        return
    end
end

is_log = getappdata(fig, 'imgislog');

m_lin = uimenu(f, 'Label', 'Linear scale', 'Tag', 'imgMenuLinear', 'Callback', @set_linear);
m_log = uimenu(f, 'Label', 'Log scale',    'Tag', 'imgMenuLog',    'Callback', @set_log);
if is_log
    m_lin.Checked = 'off';
    m_log.Checked = 'on';
else
    m_lin.Checked = 'on';
    m_log.Checked = 'off';
end
uimenu(f, 'separator', 'on', 'Label', 'Set color limits...', 'Callback', @set_clim);
uimenu(f, 'Label', 'Reset color limits', 'Callback', @reset_clim);

    function set_linear(src, ~)
        fig = ancestor(src, 'figure');
        im = findobj(fig, 'type', 'image');
        if isempty(im); return; end
        raw = getappdata(fig, 'imgraw');
        if isempty(raw); return; end
        im(1).CData = raw;
        setappdata(fig, 'imgislog', false);
        set(findobj(fig, 'tag', 'imgMenuLinear'), 'Checked', 'on');
        set(findobj(fig, 'tag', 'imgMenuLog'),    'Checked', 'off');
    end

    function set_log(src, ~)
        fig = ancestor(src, 'figure');
        im = findobj(fig, 'type', 'image');
        if isempty(im); return; end
        raw = getappdata(fig, 'imgraw');
        if isempty(raw); return; end
        im(1).CData = log10(abs(raw));
        setappdata(fig, 'imgislog', true);
        set(findobj(fig, 'tag', 'imgMenuLinear'), 'Checked', 'off');
        set(findobj(fig, 'tag', 'imgMenuLog'),    'Checked', 'on');
    end

    function set_clim(src, ~)
        fig = ancestor(src, 'figure');
        im = findobj(fig, 'type', 'image');
        if isempty(im); return; end
        ax = ancestor(im(1), 'axes');
        current = ax.CLim;
        answer = inputdlg({'Min:', 'Max:'}, 'Color limits', 1, ...
            {sprintf('%.4g', current(1)), sprintf('%.4g', current(2))});
        if isempty(answer); return; end
        new_min = str2double(answer{1});
        new_max = str2double(answer{2});
        if ~isnan(new_min) && ~isnan(new_max) && new_min < new_max
            ax.CLim = [new_min, new_max];
        end
    end

    function reset_clim(src, ~)
        fig = ancestor(src, 'figure');
        im = findobj(fig, 'type', 'image');
        if isempty(im); return; end
        ax = ancestor(im(1), 'axes');
        set(ax, 'CLimMode', 'auto');
    end

end
