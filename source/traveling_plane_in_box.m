function [pts, P] = traveling_plane_in_box(varargin)
% For a given vector, it defines a plane that include position P.
% the position P is defined as a position at a distance
% dmax*fractional_position away from O along the vector.
%
% Use this function to draw a plane to intersect a box.
% see drawvoxels.m and drawPlane3d.m
%
% 1/22/2022
% Byeongdu Lee

if numel(varargin)>3
    planevector = varargin{1};
    O = varargin{2};
    fractional_position = varargin{3};
    dmax = varargin{4};

    if numel(varargin)<5
        box = [xlim, ylim, zlim];
    else
        box = varargin{5};
    end
    % now, the plane to cut
    P = O + planevector/norm(planevector)*fractional_position*dmax;
else
    planevector = varargin{1};
    P = varargin{2};
    box = varargin{3};
end
    % Real space coordinate = fractional coordinate*mat;
    % fractional coordinate = Realspace coordinate * inv(mat)';
    pts = intersect_planebox(planevector, P, box); % plane confined by the box.        

    % putting pts in correct order.....
    plane = createPlane(P, planevector);
    % the two spanning lines of the plane
    d1 = plane(:, [1:3 4:6]);
    d2 = plane(:, [1:3 7:9]);

    % position of intersection points in plane coordinates
    u1 = linePosition3d(pts, d1);
    u2 = linePosition3d(pts, d2);

    % reorder vertices in the correct order
    ind = convhull(u1, u2);
    ind = ind(1:end-1);
    pts = pts(ind, :);
