function obj2 = polyhedron_area_volume(obj)
if isfield(obj, 'vertices')
    vert_name = 'vertices';
    face_name = 'faces';
else
    vert_name = 'Vertices';
    face_name = 'Faces';
end
[~, V] = convhull(obj.(vert_name));
try
    obj2 = obj;
    obj2.volume = V;
catch
    obj2 = [];
    obj2.volume = V;
end
face = obj.(face_name);
vert = obj.(vert_name);
if iscell(face)
    Nface = numel(face);
else
    Nface = size(face, 1);
end
S = 0;
facearea = zeros(Nface, 1);
for k=1:Nface
    if iscell(face)
        plane = vert(face{k}, :);
    else
        plane = vert(face(k, :), :);
    end
    facearea(k) = polygonArea3d(plane);
end
obj2.faceArea = facearea;
obj2.area = sum(facearea);
obj2.IsoperimetricQuotient = 36*pi*V^2/obj2.area^3;