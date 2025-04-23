function h = add_2dgrid_on3d(varargin)
% add_2dgrid_on3d('x', 0.1)
% add_2dgrid_on3d('z', 0.1, 'linewidth',1, 'linestyle', ':', 'color', 'k')
xl = xlim;
yl = ylim;
zl = zlim;
switch varargin{1}
    case 'x'
        ax1 = yl;
        ax2 = zl;
        ax3 = [varargin{2},varargin{2}];
    case 'y'
        ax1 = xl;
        ax2 = zl;
        ax3 = [varargin{2},varargin{2}];
    case 'z'
        ax1 = xl;
        ax2 = yl;
        ax3 = [varargin{2},varargin{2}];
end

axl1 = fix(ax1(1)):1:fix(ax1(2));
axl2 = fix(ax2(1)):1:fix(ax2(2));

h = [];
for i=1:numel(axl2)
    switch varargin{1}
        case 'z'
            k = plot3(ax1, [axl2(i), axl2(i)], ax3);
        case 'y'
            k = plot3(ax1, ax3, [axl2(i), axl2(i)]);
        case 'x'
            k = plot3(ax3, [axl2(i), axl2(i)], ax1);
    end
    h = [h, k];
end
for i=1:numel(axl1)
    switch varargin{1}
        case 'z'
            k = plot3([axl1(i), axl1(i)], ax2, ax3);
        case 'y'
            k = plot3([axl1(i), axl1(i)], ax3, ax2);
        case 'x'
            k = plot3(ax3, ax2, [axl1(i), axl1(i)]);
    end
    h = [h, k];
end
if numel(varargin)>2
    for i=3:2:numel(varargin)-1
        set(h, varargin{i}, varargin{i+1});
    end
end