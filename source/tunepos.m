function cenms = tunepos(varargin)

%v = [x(1), y(1), z(1)];
if numel(varargin) >= 2
    vox = varargin{1};
    pos = varargin{2};
end

th = 0;

if numel(varargin) == 3
    th = varargin{3};
end
sz = size(vox);
x = pos(:,1); y = pos(:,2); z = pos(:,3);
if numel(varargin) < 5
    sz = size(vox);
    [X, Y, Z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
end
if numel(varargin) == 5
    X = varargin{3};
    Y = varargin{4};
    Z = varargin{5};
end
if numel(varargin) == 6
    th = varargin{3};
    X = varargin{3};
    Y = varargin{4};
    Z = varargin{5};
end
ds = 2;
cenms = zeros(size(pos));
%cenms2 = cenms;
%value = [];
totalN = prod(sz);
isstarted = false;
for i=1:size(pos, 1)
    if ~isstarted
        tic
    end

    [v, xs, ys, zs] = extract_sub(vox, x(i), y(i), z(i), X, Y, Z, sz, ds);

%     % method 1:
%     %cenms(i, :) = centerofmass(v, xs, ys, zs);
%     %cen_r = round(cen);
%     % method 2: % this method works well....
     [~, ind] = max(v(:));
     cen_r = [xs(ind), ys(ind), zs(ind)];
%      cenms2(i, :) = cen_r;
     [v, xs, ys, zs, xr, yr, zr] = extract_sub(vox, cen_r(1), cen_r(2), cen_r(3), X, Y, Z, sz, ds);
     [mv, ind] = max(v(:));
     if mv < th
         cen_r = [0, 0,0];
     else
         cen_r = [xs(ind), ys(ind), zs(ind)];
         if cen_r(1)==xr(1) | cen_r(1)==xr(end) | cen_r(2)==yr(1) | cen_r(2)==yr(end) |cen_r(3)==zr(1) | cen_r(3)==zr(end)
             cen_r = [0, 0,0];
         end
     end
%     % method 3: fitting.
%      cenms2(i, :) = cp;
%     [cen_r, val, sig] = fit_position_in_voxel(v, xs, ys, zs, ...
%         cen_r(1), cen_r(2), cen_r(3));
% %     pos(i, :)
% %     cp
% %     cen_r
% %     sqrt(sum((cen_r-cp).^2)), sig, (mv-val)
%     if mod(i, fix(size(pos, 1)/10)) == 0
%         fprintf('Fit done. %i/%i\n', i, size(pos, 1));
%     end
    cenms(i, :) = cen_r;
    if ~isstarted
        t = toc;
    end
    if mod((i-1)/totalN*100, 10) == 0
        remtime = (totalN - i)*t;
        fprintf('Remaining time %0.3f seconds\n', remtime)
    end
end
t = cenms(:,1)==0 &cenms(:,2)==0 &cenms(:,3)==0;
cenms(t, :) = [];

function [v, xs, ys, zs, xr, yr, zr] = extract_sub(vox, x, y, z, X, Y, Z, sz, ds)
    if x>ds
        sxr = x-ds;
    else
        sxr = x;
    end
    if x<sz(1)-ds+1
        exr = x+ds;
    else
        exr = x;
    end
    if y>ds
        syr = y-ds;
    else
        syr = y;
    end
    if y<sz(2)-ds+1
        eyr = y+ds;
    else
        eyr = y;
    end
    if z>ds
        szr = z-ds;
    else
        szr = z;
    end
    if z<sz(3)-ds+1
        ezr = z+ds;
    else
        ezr = z;
    end
    xr = sxr:exr;
    yr = syr:eyr;
    zr = szr:ezr;
    
    xs = X(xr, yr, zr);
    ys = Y(xr, yr, zr);
    zs = Z(xr, yr, zr);
    v = vox(xr, yr, zr);
