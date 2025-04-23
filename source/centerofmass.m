function center = centerofmass(vox, x, y, z)
% Calculate the center of mass of an array.
% center = centerofmass(vox)
% center = centerofmass(vox, x, y, z)
% center = centerofmass(vox, x, y)
% output:
%   when vox is 3D,
%   center should be [xc, yc, zc]
%   , where xc is the center of the 1st dimension
%   , where yc is the center of the 1st dimension
%   , where zc is the center of the 1st dimension
% % for example, if
% b(:,:,1) =
%           10           1           6
%        10000           5           7
%            4           9           2
% b(:,:,2) =
%      8     1     6
%      3     5     7
%      4     9     2
% b(:,:,3) =
%      8     1     6
%      3     5     7
%      4     9     2
% centerofmass(b)
% ans = 
% 1.9998    1.0133    1.0133
%

if nargin<2
    sz = size(vox);
    if ndims(vox)==3
        [x, y, z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
    else
        [x, y] = ndgrid(1:sz(1), 1:sz(2));
    end
end
voxsum = sum(vox, 'all');
center(1) = sum(vox.*x, 'all')/voxsum;
center(2) = sum(vox.*y, 'all')/voxsum;
if ndims(vox)==3
    center(3) = sum(vox.*z, 'all')/voxsum;
end