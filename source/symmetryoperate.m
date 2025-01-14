function pos = symmetryoperate(sginfo, frac_coordinate, varargin)
% frac_coordinate = symmetryoperate(sginfo, frac_coordinate)
% hkl = symmetryoperate(sginfo, hkl, 'hkl', true)
% only unique pairs will be returned.

if ~isstruct(sginfo)
    if ischar(sginfo)
        sginfo = sgroup(sginfo);
    end
end
sg = sginfo;
ishkl = false;
if numel(varargin)>=2
    for k=1:2:numel(varargin)
        switch varargin{k}
            case 'hkl'
                ishkl = varargin{k+1};
        end
    end
end
SymM = sg.SymMatrices(:, 1:4, 1:sg.NoSymMatrices);

if size(sg.SymMatrices, 3)<24
    NoLatticeCenteringVector = 1; % R centering;
else
    NoLatticeCenteringVector = sg.NoLatticeCenteringVector;
end

% cm is the full list of particles
p_pos = frac_coordinate;
pos = [];
% first SymM is eye, so skip
for m = 2:size(SymM, 3)
    vt = [];
    if ~ishkl
        v2 = p_pos*SymM(:,1:3,m)' + SymM(:,4,m)'/12;
        if sg.Flag_1bar == 1
            v2 = [v2;-p_pos*SymM(:,1:3,m)' + SymM(:,4,m)'/12];
        end
        for lc = 1:NoLatticeCenteringVector
            vt = [vt;v2(:, 1:3) + sg.LatticeCenteringVector(1:3,lc)'];
        end
        t = vt < 0;
        vt(t) = vt(t)+1;
        t = vt >= 1;
        vt(t) = vt(t)-1;
    else
        v2 = p_pos*SymM(:,1:3,m)';
        if sg.Flag_1bar == 1
            v2 = [v2;-v2];
        end
        vt = v2;
    end
    pos = [pos;vt];
end
pos = unique_m(pos);

%pos = fillunitcell(pos, 1);
