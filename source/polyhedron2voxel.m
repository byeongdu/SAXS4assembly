function [vox, X, Y, Z] = polyhedron2voxel(obj, varargin)
resolution = 1;
% xmin = 0;
% xmax = 0;
% ymin = 0;
% ymax = 0;
% zmin = 0;
% zmax = 0;
shellradius = 0;
fillsquaretile = false;
filloctagontile = false;
if numel(varargin)>0
    for i=1:2:numel(varargin)
        eval(sprintf('%s = %f;', varargin{i}, varargin{i+1}));
    end
    if shellradius ~= 0 
        fprintf("shell radius = %0.3f\n", shellradius)
    end
    v = min(obj.vertices-shellradius);
    xmin = v(1); ymin = v(2); zmin = v(3);
    v = max(obj.vertices+shellradius);
    xmax = v(1); ymax = v(2); zmax = v(3);

else
    v = min(obj.vertices);
    xmin = v(1); ymin = v(2); zmin = v(3);
    v = max(obj.vertices);
    xmax = v(1); ymax = v(2); zmax = v(3);
end

xmin = xmin-resolution;
ymin = ymin-resolution;
zmin = zmin-resolution;
xmax = xmax+resolution;
ymax = ymax+resolution;
zmax = zmax+resolution;

x = xmin:resolution:xmax;
y = ymin:resolution:ymax;
z = zmin:resolution:zmax;
[X, Y, Z] = ndgrid(x, y, z);
edge = [];
tile = [];
if iscell(obj.faces)
    nf = numel(obj.faces);
else
    nf = size(obj.faces, 1);
end
for i=1:nf
    if iscell(obj.faces)
        f = obj.faces{i};
    else
        f = obj.faces(i, :);
    end
    for j=1:numel(f)
        v1 = obj.vertices(f(j), :);
        if j == numel(f)
            v2 = obj.vertices(f(1), :);
        else
            v2 = obj.vertices(f(j+1), :);
        end
        edge = [edge;v1, v2];
    end
end
sz = size(X);
vox = zeros(sz);
X = X(:);
dist = zeros(size(X));
for i=1:numel(X)
    dist(i) = min(distancePointEdge3d([X(i), Y(i), Z(i)], edge));
    if dist(i) <= shellradius
        vox(i) = 1;
    end
end
if fillsquaretile || filloctagontile
    inP = zeros(size(X));
    for i=1:nf
        if iscell(obj.faces)
            f = obj.faces{i};
        else
            f = obj.faces(i, :);
        end
        v = obj.vertices(f, :);
        if (numel(f) ==4) && (fillsquaretile)
            nv = faceNormal(v, [1,2,3,4]);
            vface = [v-nv*shellradius;v+nv*shellradius];
            faces = [1,2,3,4;
                5,6,7,8;
                2,6,7,3;
                1,5,6,2;
                4,8,5,1;
                7,3,4,8;
                6,2,3,7];
            sqfaceobj.vertices = vface;
            sqfaceobj.faces = faces;
            in = inPolyhedron2(sqfaceobj, [X(:), Y(:), Z(:)]);
            inP = inP | in;
        end
        if (numel(f) ==8) && filloctagontile
            nv = faceNormal(v, [1,2,3]);
            vface = [v-nv*shellradius;v+nv*shellradius];
            faces = findfaces(vface);
%             faces = [1,2,3,4;
%                 5,6,7,8;
%                 2,6,7,3;
%                 1,5,6,2;
%                 4,8,5,1;
%                 7,3,4,8;
%                 6,2,3,7];
            sqfaceobj.vertices = vface;
            sqfaceobj.faces = faces;
            in = inPolyhedron2(sqfaceobj, [X(:), Y(:), Z(:)]);
            inP = inP | in;
        end
    end
    vox(inP) = 1;
end

X = reshape(X, sz);
Y = reshape(Y, sz);
Z = reshape(Z, sz);
