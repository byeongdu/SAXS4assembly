function projection_result = proj4vec(voxel_map, projection_direction)
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
rotation_axis = cross([0, 0, 1], projection_direction);
rotation_axis = rotation_axis / norm(rotation_axis);

% Calculate the angle of rotation
theta = acos(dot([0, 0, 1], projection_direction));

% Create the rotation matrix using the axis-angle formula (Rodrigues' rotation formula)
K = [0, -rotation_axis(3), rotation_axis(2); 
     rotation_axis(3), 0, -rotation_axis(1); 
     -rotation_axis(2), rotation_axis(1), 0];

R = eye(3) + sin(theta)*K + (1-cos(theta))*(rotation_axis' * rotation_axis);

% Step 3: Interpolate voxel map to rotated coordinates
        [X, Y, Z] = ndgrid(1:nx, 1:ny, 1:nz);
        original_coords = [X(:), Y(:), Z(:)]'; % Original coordinates
        rotated_coords = R * original_coords;  % Rotated coordinates
        
        % Interpolation (griddata or interpn for large arrays)
        rotated_voxel_map = interpn(X, Y, Z, voxel_map, ...
            reshape(rotated_coords(1, :), nx, ny, nz), ...
            reshape(rotated_coords(2, :), nx, ny, nz), ...
            reshape(rotated_coords(3, :), nx, ny, nz), ...
            'linear', 0); % Fill with 0 for out-of-bound coordinates
        
        % Step 4: Project along z-axis (spatial domain summation)
        projection_result = sum(rotated_voxel_map, 3);
