function [Vox, x, y, z, vox] = fq2vox(Fq, dq, varargin)
% vox = fq2vox(Fq, dq)
% Calculate rho from Fq using Fourier Transform.
% vox is 3D density, not a unit cell. And thus rectangular or cubic box.
% dq is the length of each Fq.
%
% Byeongdu Lee
% 12/26/2019

maxsize = 256;

for i=1:2:numel(varargin)-1
    switch varargin{i}
        case 'maxsize'
            maxsize = varargin{i+1};
    end
end

Na = 0;
Nb = 0;
Nc = 0;
z = [];
dimdt = ndims(Fq);
if dimdt == 3
    if size(Fq, 1) <maxsize
        if rem(size(Fq, 1), 2)>0
            Fq(end, :, :) = [];
        end
        Na = (maxsize - size(Fq, 1))/2;
    end
    if size(Fq, 3) <maxsize
        if rem(size(Fq, 3), 2)>0
            Fq(:, :, end) = [];
        end
        Nc = (maxsize - size(Fq, 3))/2;
    end
    if size(Fq, 2) <maxsize
        if rem(size(Fq, 2), 2)>0
            Fq(:, end, :) = [];
        end
        Nb = (maxsize - size(Fq, 2))/2;
    end
else
    if size(Fq, 1) <maxsize
        if rem(size(Fq, 1), 2)>0
            Fq(end, :) = [];
        end
        Na = (maxsize - size(Fq, 1))/2;
    end
    if size(Fq, 2) <maxsize
        if rem(size(Fq, 2), 2)>0
            Fq(:, end) = [];
        end
        Nb = (maxsize - size(Fq, 2))/2;
    end
end
% not sure whether zero padding is necessary.
if dimdt == 3
    Fq = padarray(Fq, [Na, Nb, Nc], 0, 'both');
else
    Fq = padarray(Fq, [Na, Nb], 0, 'both');
end

vox = ifftn(Fq);
vox = ifftshift(vox, 1);
vox = ifftshift(vox, 2);
if dimdt == 3
    vox = ifftshift(vox, 3);
end
%Iq = abs(F).^2;
Vox = vox;
Vox = Vox./max(Vox, [], 'all');
if nargin > 1
    if (dimdt == 3) && (numel(dq) == 1)
        dq(2) = dq(1);
        dq(3) = dq(1);
    end
    if (dimdt == 2) && (numel(dq) == 1)
        dq(2) = dq(1);
    end
    x = 1:size(Vox, 2);
    x = x - size(Vox, 2)/2-1;
    x = x*2*pi/dq(1)/size(Vox, 2);
    y = 1:size(Vox, 1);
    y = y - size(Vox, 1)/2-1;
    y = y*2*pi/dq(2)/size(Vox, 1);
    if dimdt == 3
        z = 1:size(Vox, 3);
        z = z - size(Vox, 3)/2-1;
        z = z*2*pi/dq(3)/size(Vox, 3);
    end
end
%cnt = size(Iq)/2+1;
%dd = Iq(cnt(1)-20:cnt(1)+20,cnt(2)-20:cnt(2)+20,cnt(3)-20:cnt(3)+20);