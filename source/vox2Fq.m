function FF = vox2Fq(Vox, qx, qy, qz, voxel_size)
% FF = vox2Fq(Vox, qx, qy, qz, voxel_size)
% dimension of qx, qy and qz should be the same.
% Calculate the scattering amplitude from a  Vox.
%
% Byeongdu Lee
% 1/2/2024

siz = size(Vox);
nd = ndims(Vox);
if nd == 1
    x = linspace(1, siz(1), siz(1));
    y = 0;
    z = 0;
end
if nd == 2
    x = linspace(1, siz(2), siz(2));
    y = linspace(1, siz(1), siz(1));
    z = 0;
end
if nd == 3
    x = linspace(1, siz(2), siz(2));
    y = linspace(1, siz(1), siz(1));
    z = linspace(1, siz(3), siz(3));
end
[Xax, Yax, Zax] = meshgrid(x, y, z);

if nargin < 5
    voxel_size = 1;
end
Xax = Xax*voxel_size;
Yax = Yax*voxel_size;
Zax = Zax*voxel_size;
FF = zeros(size(qx));
t = Vox==0;
Vox(t) = [];
Xax(t) = [];
Yax(t) = [];
Zax(t) = [];
for i=1:numel(qx)
    F= Vox.*exp(-sqrt(-1)*(qx(i)*Xax+qy(i)*Yax+qz(i)*Zax));
    FF(i) = sum(F(:));
end
FF = FF/numel(Vox);