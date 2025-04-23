function [intensity, x, y, z] = vox_hkl_integrate(vox, pos, X, Y, Z, ds)
sz = size(vox);
x = pos(:,1); y = pos(:,2); z = pos(:,3);
intensity = zeros(size(pos, 1), 1);
qv = [X(:), Y(:), Z(:)];
%value = [];
totalN = prod(sz);
isstarted = false;
for i=1:size(pos, 1)
    d = sqrt(sum((qv-pos(i, :)).^2, 2));
    [~, ind] = min(d);
    [xind, yind, zind] = ind2sub(sz, ind);
    x(i) = xind;
    y(i) = yind;
    z(i) = zind;
    if ~isstarted
        tic
    end
    if x(i)>ds
        sxr = x(i)-ds;
    else
        sxr = x(i);
    end
    if x(i)<sz(1)-ds+1
        exr = x(i)+ds;
    else
        exr = x(i);
    end
    if y(i)>ds
        syr = y(i)-ds;
    else
        syr = y(i);
    end
    if y(i)<sz(2)-ds+1
        eyr = y(i)+ds;
    else
        eyr = y(i);
    end
    if z(i)>ds
        szr = z(i)-ds;
    else
        szr = z(i);
    end
    if z(i)<sz(3)-ds+1
        ezr = z(i)+ds;
    else
        ezr = z(i);
    end
    xr = sxr:exr;
    yr = syr:eyr;
    zr = szr:ezr;
    
    xs = X(xr, yr, zr);
    ys = Y(xr, yr, zr);
    zs = Z(xr, yr, zr);
    
    % read the sub-cube and fill Nan with average value.
    v = vox(xr, yr, zr);
    if sum(isnan(v)) == numel(v(:))
        intensity(i) = nan;
    elseif sum(isnan(v)) == 0
        intensity(i) = sum(v(2:end-1, 2:end-1, 2:end-1), 'all');
    else
        f = 1;
        while f>0
            [v, f] = fillNaN(v);
        end
        intensity(i) = sum(v(2:end-1, 2:end-1, 2:end-1), 'all');
    end
    center = centerofmass(v(2:end-1, 2:end-1, 2:end-1), ...
        xs(2:end-1, 2:end-1, 2:end-1), ...
        ys(2:end-1, 2:end-1, 2:end-1), ...
        zs(2:end-1, 2:end-1, 2:end-1));
    x(i) = center(1);
    y(i) = center(2);
    z(i) = center(3);
    if ~isstarted
        t = toc;
    end
    if mod((i-1)/totalN*100, 10) == 0
        remtime = (totalN - i)*t;
        fprintf('Remaining time %0.3f seconds\n', remtime)
    end
end
