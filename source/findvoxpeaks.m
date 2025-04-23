function [cenms, val] = findvoxpeaks(varargin)
% c = findvoxpeaks(img, ds)
% c = findvoxpeaks(img, X, Y, Z, ds)
% output:
%    c = [position_dim1, position_dim2, position_dim3]
%        position_dim1 : row (therefore y)
%        position_dim2 : column (therefore x)
%        position_dim3 : z
% To read values at the position,
%   ind = sub2ind(size(img), c(:,1), c(:,2), c(:,3));
%   val = img(ind);
% To plot the position:
%   plot3(c(:,2), c(:,1), c(:,3), 'ro')
%

%[~, ind] = max(img(:));
%[aa, bb, cc] = ind2sub(size(dt), ind);
if numel(varargin) < 5
    ds = 2;
end
img = varargin{1};
sz = size(img);
if numel(varargin)==2
    [x, y, z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
    ds = varargin{2};
else
    x = varargin{2};
    y = varargin{3};
    z = varargin{4};
    ds = varargin{5};
end
cenms = [];
value = [];
totalN = prod(sz);
isstarted = false;
for i=1:size(img, 1)
    for j=1:size(img, 2)
        for k=1:size(img, 3)
            if ~isstarted
                tic
            end
            if i>ds
                sxr = i-ds;
            else
                sxr = i;
            end
            if i<sz(1)-ds+1
                exr = i+ds;
            else
                exr = i;
            end
            if j>ds
                syr = j-ds;
            else
                syr = j;
            end
            if j<sz(2)-ds+1
                eyr = j+ds;
            else
                eyr = j;
            end
            if k>ds
                szr = k-ds;
            else
                szr = k;
            end
            if k<sz(3)-ds+1
                ezr = k+ds;
            else
                ezr = k;
            end
            xr = sxr:exr;
            yr = syr:eyr;
            zr = szr:ezr;
            xs = x(xr, yr, zr);
            ys = y(xr, yr, zr);
            zs = z(xr, yr, zr);
            v = img(xr, yr, zr);
            % method 1:
            %cen = centerofmass(v, xs, ys, zs);
            %cen_r = round(cen);
            % method 2:
            %[val, ind] = max(v(:));
            %cen_r = [xs(ind), ys(ind), zs(ind)];
            % method 3: fitting.
            [~, ind] = max(v(:));
            x0 = xs(ind);
            y0 = ys(ind);
            z0 = zs(ind);
            [cen_r, val] = fit_position_in_voxel(v, xs, ys, zs, x0, y0, z0);
%            if cen_r(1) == i & cen_r(2) == j & cen_r(3) ==k
                    cenms = [cenms;cen_r];
                    value = [value; val];
%            end
            if ~isstarted
                t = toc;
            end
            if mod((i*j*k-1)/totalN*100, 10) == 0
                remtime = (totalN - i*j*k)*t;
                fprintf('Remaining time %0.3f seconds\n', remtime)
            end
        end
    end
end

% 3d ball peak