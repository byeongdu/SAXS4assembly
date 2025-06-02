function [Iq, qx, qy, qz, F, F2] = vox2Iq(dt, vox_size, varargin)
% [Iq, qx, qy, qz, F] = vox2Iq(dt, vox_size)
% qx, qy, qz are ndgrid.
% Reciprcal space of dt using Fourier Transform.
% Iq = vox2Iq(dt)
% dt is 3D density, not a unit cell. And thus rectangular or cubic box.
% vox_size is the length of each voxel that is a cubic.
%
% Byeongdu Lee
% 12/26/2019

maxsize = 256;
F2 = [];
isRetF2 = false;
if nargout==6
    isRetF2 = true;
end
backg = 0;
for i=1:2:numel(varargin)-1
    switch varargin{i}
        case 'maxsize'
            maxsize = varargin{i+1};
        case 'background'
            if ischar(varargin{i+1})
                if strcmp(varargin{i+1}, 'mean')
                    backg = nanmean(dt(:));
                end
                if strcmp(varargin{i+1}, 'min')
                    backg = min(dt(:));
                end
            else
                backg = varargin{i+1};
            end
    end
end

Na = 0;
Nb = 0;
Nc = 0;
qz = [];

dimdt = ndims(dt);

if dimdt == 3
    if size(dt, 1) <maxsize
        if rem(size(dt, 1), 2)>0
            dt(end, :, :) = [];
        end
        Na = (maxsize - size(dt, 1))/2;
    end
    if size(dt, 3) <maxsize
        if rem(size(dt, 3), 2)>0
            dt(:, :, end) = [];
        end
        Nc = (maxsize - size(dt, 3))/2;
    end
    if size(dt, 2) <maxsize
        if rem(size(dt, 2), 2)>0
            dt(:, end, :) = [];
        end
        Nb = (maxsize - size(dt, 2))/2;
    end
else
    if size(dt, 1) <maxsize
        if rem(size(dt, 1), 2)>0
            dt(end, :) = [];
        end
        Na = (maxsize - size(dt, 1))/2;
    end
    if size(dt,2)>1
        if size(dt, 2) <maxsize
            if rem(size(dt, 2), 2)>0
                dt(:, end) = [];
            end
            Nb = (maxsize - size(dt, 2))/2;
        end
    end
end
%dtn = dt-mean(dt(:));
dtn = dt-backg;
if dimdt == 3
    dtn = padarray(dtn, [Na, Nb, Nc], 0, 'both');
    if isRetF2
        dt = padarray(dt, [Na, Nb, Nc], 0, 'both');
    end
else
    dtn = padarray(dtn, [Na, Nb], 0, 'both');
    if isRetF2
        dt = padarray(dt, [Na, Nb], 0, 'both');
    end
end

if isRetF2
    F2 = fftn(dt);
    F2 = fftshift(F2, 1);
    F2 = fftshift(F2, 2);
    if dimdt == 3
        F2 = fftshift(F2, 3);
    end
end

F = fftn(dtn);
F = fftshift(F, 1);
F = fftshift(F, 2);
if dimdt == 3
    F = fftshift(F, 3);
end
Iq = abs(F).^2;
Iq = Iq./max(Iq, [], 'all');
if nargin > 1
    if (dimdt == 3) && (isscalar(vox_size))
        vox_size(2) = vox_size(1);
        vox_size(3) = vox_size(1);
    end
    if (dimdt == 2) && (isscalar(vox_size))
        vox_size(2) = vox_size(1);
    end
    qx = 1:size(Iq, 1);
    qx = qx - size(Iq, 1)/2-1;
    qx = qx*2*pi/vox_size(1)/size(Iq, 1);
    qy = 1:size(Iq, 2);
    qy = qy - size(Iq, 2)/2-1;
    qy = qy*2*pi/vox_size(2)/size(Iq, 2);
    if dimdt == 3
        qz = 1:size(Iq, 3);
        qz = qz - size(Iq, 3)/2-1;
        qz = qz*2*pi/vox_size(3)/size(Iq, 3);
    end
end
%cnt = size(Iq)/2+1;
%dd = Iq(cnt(1)-20:cnt(1)+20,cnt(2)-20:cnt(2)+20,cnt(3)-20:cnt(3)+20);