function map = symmetryoperate2voxel(map, sginfo, frac_start)
% frac_start : when voxel's fractional coordinate ranges from 0 to 1,
% frac_start = 0
%   for example, when use fft.
% if voxel does not symmetric (0 to 1), use option 1.
if nargin<3 % voxel's fractional coordinate ranges from 0 to 1.
    frac_start = 1;
end

SymM = sginfo.SymMatrices(:, 1:4, 1:sginfo.NoSymMatrices);

if size(sginfo.SymMatrices, 3)<24
    NoLatticeCenteringVector = 1; % R centering;
else
    NoLatticeCenteringVector = sginfo.NoLatticeCenteringVector;
end

    siz = size(map);
    [Xp, Yp, Zp] = ndgrid(1:siz(1), 1:siz(2), 1:siz(3));
    if frac_start == 0
        frac_coord = [(Xp(:)-1)/(siz(1)-1), (Yp(:)-1)/(siz(2)-1), (Zp(:)-1)/(siz(3)-1)];
    else
        frac_coord = [Xp(:)/siz(1), Yp(:)/siz(2), Zp(:)/siz(3)];
    end
    for lc = 1:NoLatticeCenteringVector
        lv = sginfo.LatticeCenteringVector(1:3,lc)';
        frac_coord2 = sym_centering(frac_coord, lv);
        map = apply_sym(map, siz, frac_coord2, Xp, Yp, Zp, frac_start);
    end
    if sginfo.Flag_1bar == 1
        frac_coord2 = sym_inv(frac_coord);
        map = apply_sym(map, siz, frac_coord2, Xp, Yp, Zp, frac_start);
    end
    for m = 2:size(SymM, 3)
        frac_coord2 = sym_op(SymM, frac_coord, m);
        map = apply_sym(map, siz, frac_coord2, Xp, Yp, Zp, frac_start);
    end
    
    %frac_coord2 = symmetryoperate(sginfo, frac_coord);
    % frac_coordinate to array indices

end

function map = apply_sym(map, siz, frac_coord2, Xp, Yp, Zp, frac_start)
    X = frac_coord2(:,1);
    Y = frac_coord2(:,2);
    Z = frac_coord2(:,3);
    if frac_start == 0 % ex) 0, 0.1, 0.2, ..., 1  (11 voxels)
        X = round(X*(siz(1)-1)+1);
    else % ex) 0.1, 0.2, ...., 1 (10 voxels)
        X(X==0)=1;
        X = round(X*siz(1));
    end
    if frac_start == 0 % ex) 0, 0.1, 0.2, ..., 1  (11 voxels)
        Y = round(Y*(siz(2)-1)+1);
    else % ex) 0.1, 0.2, ...., 1 (10 voxels)
        Y(Y==0)=1;
        Y = round(Y*siz(1));
    end
    if frac_start == 0 % ex) 0, 0.1, 0.2, ..., 1  (11 voxels)
        Z = round(Z*(siz(3)-1)+1);
    else % ex) 0.1, 0.2, ...., 1 (10 voxels)
        Z(Z==0)=1;
        Z = round(Z*siz(1));
    end
    map0 = map;
    ind = sub2ind(siz, Xp(:), Yp(:), Zp(:));
    ind2 = sub2ind(siz, X(:), Y(:), Z(:));
    map = (map0(ind)+map0(ind2))/2;
    map = reshape(map, siz);
end
function frac = sym_op(SymM, frac, m)
    frac = frac*SymM(:,1:3,m)' + SymM(:,4,m)'/12;
    t = frac < 0;
    frac(t) = frac(t)+1;
    t = frac >= 1;
    frac(t) = frac(t)-1;
end
function frac = sym_inv(frac)
    frac = -frac;
    t = frac < 0;
    frac(t) = frac(t)+1;
    t = frac >= 1;
    frac(t) = frac(t)-1;
end
function frac = sym_centering(frac, LatticeCenteringVector)
    frac = frac + LatticeCenteringVector;
    t = frac < 0;
    frac(t) = frac(t)+1;
    t = frac >= 1;
    frac(t) = frac(t)-1;
end
