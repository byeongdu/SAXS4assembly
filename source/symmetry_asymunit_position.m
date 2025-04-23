function Ind_out = symmetry_asymunit_position(map, sginfo)
%function map = symmetry_asymunit_position(X, Y, Z, Nvoxel, sginfo)
% Note : fractional coordinate should be > 0 and <= 1
% ex)   X = linspace(0, 1, Nvoxel+1); X(1) = [];
%       Xi = X*Nvoxel;
% map should be cubic volume with the number of voxel is Nvoxel along X.
%   for example, when use fft.
% if voxel does not symmetric (0 to 1), use option 1.
Nvoxel = length(map);
x = linspace(0, 1, Nvoxel+1); x(1) = [];
[X, Y, Z] = ndgrid(x,x,x);

SymM = sginfo.SymMatrices(:, 1:4, 1:sginfo.NoSymMatrices);

if size(sginfo.SymMatrices, 3)<24
    NoLatticeCenteringVector = 1; % R centering;
else
    NoLatticeCenteringVector = sginfo.NoLatticeCenteringVector;
end

frac_coord = [X(:), Y(:), Z(:)];
Ind_in = 1:numel(X(:));
Ind_out = zeros(size(Ind_in));
for ind_progress = 1:numel(X(:))
    if Ind_out(ind_progress) > 0
        continue
    end
    frac_in = frac_coord(ind_progress, :);
    for lc = 2:NoLatticeCenteringVector
        lv = sginfo.LatticeCenteringVector(1:3,lc)';
        frac = sym_centering(frac_in, lv);
        frac_in = [frac_in;frac];
        frac_in = unique_m(frac_in);
    end
    if sginfo.Flag_1bar == 1
        frac = sym_inv(frac_in);
        frac_in = [frac_in;frac];
        frac_in = unique_m(frac_in);
    end
    for m = 2:size(SymM, 3)
        frac = sym_op(SymM, frac_in, m);
        frac_in = [frac_in;frac];
        frac_in = unique_m(frac_in);
    end
    ind2 = apply_sym(Nvoxel, frac_in(2:end,:));
    Ind_out(ind2) = ind_progress;
end
    %frac_coord2 = symmetryoperate(sginfo, frac_coord);
    % frac_coordinate to array indices

end

function ind2 = apply_sym(Nvoxel, frac_coord2)
    X = frac_coord2(:,1);
    Y = frac_coord2(:,2);
    Z = frac_coord2(:,3);
    X = round(X*Nvoxel);
    Y = round(Y*Nvoxel);
    Z = round(Z*Nvoxel);
    X(X==0)=Nvoxel;
    Y(Y==0)=Nvoxel;
    Z(Z==0)=Nvoxel;
    
    %ind = sub2ind(siz, Xp(:), Yp(:), Zp(:));
    try
    ind2 = sub2ind([Nvoxel, Nvoxel, Nvoxel], X(:), Y(:), Z(:));
    catch
        [X(:), Y(:), Z(:)]
    end
    
end
function frac = sym_op(SymM, frac, m)
    frac = frac*SymM(:,1:3,m)' + SymM(:,4,m)'/12;
    t = frac <= 0;
    frac(t) = frac(t)+1;
    t = frac > 1;
    frac(t) = frac(t)-1;
    frac = frac(:, 1:3);
end
function frac = sym_inv(frac)
    frac = -frac;
    t = frac <= 0;
    frac(t) = frac(t)+1;
    t = frac > 1;
    frac(t) = frac(t)-1;
end
function frac = sym_centering(frac, LatticeCenteringVector)
    frac = frac + LatticeCenteringVector;
    t = frac <= 0;
    frac(t) = frac(t)+1;
    t = frac > 1;
    frac(t) = frac(t)-1;
end
