function vox = remove_intensity_from_vox(vox, pos, X, Y, Z, ds)
sz = size(vox);
x = pos(:,1); y = pos(:,2); z = pos(:,3);
qv = [X(:), Y(:), Z(:)];
%value = [];
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
    
    % read the sub-cube and fill Nan with average value.
    vox(xr, yr, zr)=nan;
end