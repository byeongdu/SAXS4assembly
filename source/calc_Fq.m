function FF = calc_Fq(Vox, h, k, l, cellinfo)
% Vox data : ranges from frac_x = 0 to (N-1) / N;
% FF = calc_Fq(Vox, h, k, l)
% FF = calc_Fq({frac_X, frac_Y, frac_Z, Vox}, h, k, l)
% FF = calc_Fq(Vox, qx, qy, qz, cellinfo)
% FF = calc_Fq({frac_X, frac_Y, frac_Z, Vox}, qx, qy, qz, cellinfo)
% Calculate the form factor from a 3D unit cell, Vox.
% h, k, l are miller indices
% Instead, qx, qy, qz can be given. Then, one must need cellinfo to convert
% the q values into hkl.
% Fractional coordinates for Vox can be provided. Their matrix dimension
% should be the same with that of Vox.
%
% Fhkl = Vcell/(NxNyNz) ∑_p=0 ^(Nx-1) ∑_q=0 ^(Ny-1) ∑_r=0 ^(Nz-1) ρ(xp, yq,
% zr)*exp (+2πi (hxp +kyq +lzr)))
%
% Byeongdu Lee
% 1/2/2020
Vcell = 1;
if iscell(Vox)
    frac_X = Vox{1};
    frac_Y = Vox{2};
    frac_Z = Vox{3};
    Vox = Vox{4};
else
    siz = size(Vox);
    nd = ndims(Vox);
    if nd == 1
        x = linspace(0, 1, siz(1)+1);x(end) = [];
        y = 0;
        z = 0;
    end
    if nd == 2
        x = linspace(0, 1, siz(1)+1);x(end) = [];
        y = linspace(0, 1, siz(2)+1);y(end) = [];
        z = 0;
    end
    if nd == 3
        x = linspace(0, 1, siz(1)+1);x(end) = [];
        y = linspace(0, 1, siz(2)+1);y(end) = [];
        z = linspace(0, 1, siz(3)+1);z(end) = [];
    end
    %[Xax, Yax, Zax] = meshgrid(x, y, z);
    [Xax, Yax, Zax] = ndgrid(x, y, z);
    frac_Y = reshape(Yax(:)/size(Xax, 1), size(Xax));
    frac_X = reshape(Xax(:)/size(Xax, 2), size(Xax));
    frac_Z = reshape(Zax(:)/size(Xax, 3), size(Xax));
end
if nargin == 5
    % h, k, l must be qx, qy, qz.
    hkl = [h(:), k(:), l(:)]*inv(cellinfo.recimat)';
    h = hkl(:, 1);
    k = hkl(:, 2);
    l = hkl(:, 3);
    Vcell = cellinfo.Vol;
else
    h = h(:);
    k = k(:);
    l = l(:);
end
FF = zeros(size(h));
for i=1:length(h)
    F= Vox.*exp(-2*pi*sqrt(-1)*(h(i)*frac_X+k(i)*frac_Y+l(i)*frac_Z));
    FF(i) = sum(F(:));
end
FF = FF/numel(Vox)*Vcell;