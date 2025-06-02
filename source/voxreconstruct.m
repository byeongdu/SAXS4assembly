function [img, h, U] = voxreconstruct(x, y, z, Zdir, Ydir, Vox, cellinfo, isovalue, type_of_isovalue)
% voxreconstruct(x, y, z, Zdir, Ydir, map3D, cellinfo, [isovalue], [type_of_isovalue])
% 
% When a unit cell vox is given, this program recontruct the map into new
% box that is specified by x, y and z (vectors in Cartesian coordinate) and
% Zdir and Ydir that are two hkl plane normal to [0 0 1] and [0 1 0]
% direction in the laboratory coordinate.
% <Zdir> will be perfectly lined up with Z axis of lab coordinate.
% However, <Ydir> does not necessarily, if structure is not a cubic. 
% <Ydir>'s projection along Z axis will line up with the Y axis...
% 
% isovalue can be an absolute number, then type_of_isovalue = 'absolute'
% isovalue can be a fraction, 0 to 1, then type_of_isovalue = 'fraction'
%   This is a scaled density value.
% isovalue can be 0 to 1 when type_of_isovalue = 'volume_fraction_high'
%   for top x% volume fraction.
% isovalue can be 0 to 1 when type_of_isovalue = 'volume_fraction_low'
%   for bottom x% volume fraction.
% if type_of_isovalue is not given, default will be fraction.
%
% This can be used to reconstructure a realspace film structure from a
% unit cell density map.
% map3D can be a cell containing {Xax, Yax, Zax, map3D} or just a density
% map of a unit cell, map3D.
%
% How to use:
% For example:
% x = 0:20:3500;
% y = 0:20:3500;
% z = 0:20:1000;
% Zdir = [1 2 1];Ydir = [0 1 1];
% 
% see also
% transformUC.m

[Vx, Vy, Vz, U] = slice2film(x, y, z, Zdir, Ydir);
np = cartesian2fractional([Vx(:), Vy(:), Vz(:)], cellinfo);
%np = [Vx(:), Vy(:), Vz(:)]*inv(cellinfo.mat)';
fracX = np(:,1);
fracY = np(:,2);
fracZ = np(:,3);

if iscell(Vox)
    frac_X = Vox{1};
    frac_Y = Vox{2};
    frac_Z = Vox{3};
    Vox = Vox{4};
else
    frac_X = [];
    frac_Y = [];
    frac_Z = [];
end

maxv = max(Vox, [], 'all');
minv = min(Vox, [], 'all');

if nargin<8
    isovalue = 0.2;
end
if nargin < 9
    type_of_isovalue = 'fraction';
end
switch type_of_isovalue
    case 'fraction'
        isovalue = minv + (maxv-minv)*isovalue;
        
    case 'volume_fraction_high'
        [N, edge] = histcounts(Vox(:), 1000);
        N = cumsum(N);
        N = N/N(end);N = 1-N;
        [~,indx] = min(abs(N-isovalue));
        isovalue = edge(indx);
        
    case 'volume_fraction_low'
        [N, edge] = histcounts(Vox(:), 1000);
        N = cumsum(N);
        N = N/N(end);
        [~,indx] = min(abs(N-isovalue));
        isovalue = edge(indx);
        
    case 'absolute'
        % no change
    otherwise
        isovalue = minv + (maxv-minv)*isovalue;
end

Vq = interpUC(frac_X, frac_Y, frac_Z, Vox, fracX, fracY, fracZ);
img = reshape(Vq, size(Vx));
%%
h = figure;
 [Xx, Yy, Zz] = meshgrid(x, y, z);
 if isovalue<min(Vq(:)) | isovalue>max(Vq)
     isovalue = mean(Vq);
 end
[faces, verts, cols] = isosurface(Xx, Yy, Zz, img, isovalue, Zz);
[F,V,C] = isocaps(Xx, Yy, Zz, img, isovalue);
% channels..
patch('faces',faces,...
                'vertices', verts,...
                'facevertexCData',cols,...
                'FaceColor','interp', 'edgecolor', 'none');
hold on
% Caps...
patch('faces',F,...
                'vertices', V,...
                'facevertexCData',C,...
                'FaceColor','interp', 'edgecolor', 'none');
%ylabel(sprintf('Ghkl=[%i,%i,%i] direction (nm)',Ydir), 'fontsize', 15)
%zlabel(sprintf('Ghkl=[%i,%i,%i] direction (nm)',Zdir), 'fontsize', 15)

axis image;
camlight;
lighting gouraud;
camva(8);