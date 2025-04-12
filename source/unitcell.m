function [o, box_vertices] = unitcell(varargin)
% unitcell(a1, a2, a3)
% unitcell(a1, a2, a3, O)
% unitcell(cellinfo)
% unitcell(cellinfo, O)
% a1 : vector for a axis;
% a2 : vector for b axis;
% a3 : vector for c axis;
% O : origin of the unit cell

o = [];
box_faces = [5 6 7 8; 1 2 3 4; 5 8 4 1; 6 7 3 2; 8 7 3 4; 5 6 2 1;];
box_vertices = [0     0     0
                    1     0     0
                    1     1     0
                    0     1     0
                    0     0     1
                    1     0     1
                    1     1     1
                    0     1     1
                  ];
cell_origin = [0,0,0];
if numel(varargin)<3
    if isfield(varargin{1}, 'latticevectors')
        latticevectors = varargin{1}.latticevectors;
    end
else
    latticevectors=[varargin{1};varargin{2};varargin{3}];
end
if numel(varargin)==2
    cell_origin = varargin{2};
end
if numel(varargin)==4
    cell_origin = varargin{4};
end
    
o.vertices = box_vertices;
o.vertices = o.vertices * latticevectors+cell_origin;
o.faces = box_faces;