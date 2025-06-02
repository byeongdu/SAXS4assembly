function [proj, xq, yq, zq] = proj4vec(voxel_map, projection_direction)
% Define a 3D voxel map (example)
%nx = 64; ny = 64; nz = 64;
%voxel_map = rand(nx, ny, nz) > 0.5;  % Random binary voxel map for testing
nx = size(voxel_map,1);
ny = size(voxel_map,2);
nz = size(voxel_map,3);
% Define an arbitrary projection direction as a unit vector
%projection_direction = [1, 1, 1];  % Example direction [1, 1, 1]

% Normalize the direction vector
projection_direction = projection_direction / norm(projection_direction);

% Step 1: Create a rotation matrix to align the projection direction with the z-axis
% This can be done using rotation matrices, here we use a simple approach via axis-angle rotation
% Get the axis of rotation (perpendicular to the z-axis and the projection direction)

% Calculate the angle of rotation
R = rotate_between_vectors([0,0,1], projection_direction);

% Step 3: Interpolate voxel map to rotated coordinates
        [X, Y, Z] = ndgrid(1:nx, 1:ny, 1:nz);
        center = [mean(1:nx); mean(1:ny); mean(1:nz)];
        original_coords = [X(:), Y(:), Z(:)]'; % Original coordinates
        rotated_coords = R' * (original_coords-center)+center;  % Rotated coordinates
        
        % Interpolation (griddata or interpn for large arrays)
        rotated_voxel_map = interpn(X, Y, Z, voxel_map, ...
            reshape(rotated_coords(2, :), nx, ny, nz), ...
            reshape(rotated_coords(1, :), nx, ny, nz), ...
            reshape(rotated_coords(3, :), nx, ny, nz), ...
            'linear', 0); % Fill with 0 for out-of-bound coordinates
        
        % Step 4: Project along z-axis (spatial domain summation)
        proj = mean(rotated_voxel_map, 3);
xq = reshape(rotated_coords(2, :), nx, ny, nz);
yq = reshape(rotated_coords(1, :), nx, ny, nz);
zq = reshape(rotated_coords(3, :), nx, ny, nz);        
xq = xq(:,:,1);
yq = yq(:,:,1);
zq = mean(zq,3);        
box = [1, nx, 1, ny, 1, nz];
pts = traveling_plane_in_box(projection_direction, center', box);
        cen = mean(pts);
        pts_xy = (pts-cen)*R'+cen;
%T = ([cols_f(:),rows_f(:), zeros(size(cols_f(:)))]*R'+center');
[X0, Y0] = ndgrid(1:nx, 1:ny);
t = inpolygon(X0, Y0, pts_xy(:,1), pts_xy(:,2));
proj(~t) = NaN;