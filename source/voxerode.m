function output = voxerode(varargin)
% voxerode(vox) % default fraction = 0.1
% voxerode(vox, 0.1)
% voxerode(vox, 'fraction', 0.1, 'mode', 'erosion')
% voxerode(vox, 'mode', 'erosion')
% voxerode(vox, 'mode', 'conv') % convolution mode.
% voxerode(vox, 'fraction', 0.1, 'mode', 'conv')

% IM: 2D or 3D matrix to be eroded
% r : percentage of erosion.
% 

IM = varargin{1};
mode = 'erosion';
r = 0.1;
if numel(varargin)==2
    r = varargin{2};
else
    for i=2:2:numel(varargin)
        val = varargin{i+1};
        switch varargin{i}
            case 'fraction'
                r = val;
            case 'mode'
                mode = val;
        end
    end
end

IM_size = size(IM);
rm_size = rem(IM_size,2);
pad_size = 1-rm_size;
IM = padarray(IM, pad_size, 'pre');
erode_resolution = 1;
rm = rem(size(IM)/erode_resolution,2);
%siz = size(IM);
rm = 1-rm;

if ndims(IM)==3
    SE = zeros(size(IM)+rm);
    [X,Y,Z] = ndgrid(linspace(-1,1,size(SE,1)),linspace(-1,1,size(SE,2)),linspace(-1,1,size(SE,3)));
    t = sqrt((X).^2+(Y).^2+(Z).^2)<=r;
else
    SE = zeros(size(IM)+rm);
    [X,Y] = ndgrid(linspace(-1,1,size(SE,1)),linspace(-1,1,size(SE,2)));
    t = sqrt((X).^2+(Y).^2)<=r;
end
SE(t) = 1;
output=convn(IM,SE/sum(SE(:)),'same');
if strcmp(mode, 'erosion')
    output(output<1) = 0;
end

ROIx = (1+pad_size(1)):(pad_size(1)+IM_size(1));
if ndims(IM)>1
    ROIy = (1+pad_size(2)):(pad_size(2)+IM_size(2));
else
    output = output(ROIx);
end
if ~ismatrix(IM)
    ROIz = (1+pad_size(3)):(pad_size(3)+IM_size(3));
    output = output(ROIx, ROIy, ROIz);
else
    output = output(ROIx, ROIy);
end

%output(abs(output-IM)>0) = 0;