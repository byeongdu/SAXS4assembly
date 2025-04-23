function Vq = interpUC(frac_X, frac_Y, frac_Z, map3D, newx, newy, newz)
% function Vq = interpUC(frac_X, frac_Y, frac_Z, map3D, newx, newy, newz)
% function Vq = interpUC(map3D, newx, newy, newz)
% This is a function for interplating a Unit Cell to any fractional
% coordinates.
% The fractional coordinates of map3D can be from 0 to 1.
% But, remember that 0 and 1 position should be identical.
% Normally, a UC density map does not include density at fractioanl_x = 0,
% fractonal_y = 0, and fractional_z = 0.
if isempty(frac_X)
    siz = size(map3D);
    x = linspace(0, 1, siz(2)+1);
    y = linspace(0, 1, siz(1)+1);
    z = linspace(0, 1, siz(3)+1);
    x(1) = [];
    y(1) = [];
    z(1) = [];
    [frac_X, frac_Y, frac_Z] = meshgrid(x, y, z);
%     frac_Y = Yax(:)/siz(1);
%     frac_X = Xax(:)/siz(2);
%     frac_Z = Zax(:)/siz(3);
end    

minx = floor(min(newx));
maxx = ceil(max(newx));
miny = floor(min(newy));
maxy = ceil(max(newy));
minz = floor(min(newz));
maxz = ceil(max(newz));
if minz==maxz
    minz = maxz-1;
end
if miny==maxy
    miny=maxy-1;
end
if minx==maxx
    minx = maxx-1;
end
Vq = nan*ones(size(newx));
newx = round(newx*1E4)/1E4;
newy = round(newy*1E4)/1E4;
newz = round(newz*1E4)/1E4;
for xi = minx:(maxx-1)
    for yi = miny:(maxy-1)
        for zi = minz:(maxz-1)
            t = (newx>=xi) & (newx<=xi+1);
            t = t & (newy>=yi) & (newy<=yi+1);
            t = t & (newz>=zi) & (newz<=zi+1);
            t = t & isnan(Vq);
            if sum(t) == 0
                continue;
            end
            try
            Vq(t) = interp3(frac_X+xi, frac_Y+yi, frac_Z+zi, map3D, ...
                newx(t), newy(t), newz(t)); 
            catch
            Vq(t) = interp3(frac_Y+yi, frac_X+xi, frac_Z+zi, map3D, ...
                newy(t), newx(t), newz(t)); 
            end
        end
    end
end
