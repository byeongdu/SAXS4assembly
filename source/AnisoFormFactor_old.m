function [Fqhkls, Pqhkls, Pqhkls2] = AnisoFormFactor(particles, hklsorq, cellinfo, dim)
% [Fqhkls, Pqhkls, Pqhkls2] = AnisoFormFactor(functionname, hkls, cellinfo, method)
% Output
%   Fqhkls: F(qx,qy,qz)
%   Pqhkls : sum<|F|^2>
%   Pqhkls2 : sum|<F>|^2
%
% Input
%   particles: particles, struct format
%   hklsorq: it coulbe hkls(struct), [qx, qy, qz] (N*3 matrix), or q (a colum vector)
%   cellinfo: cellinfo (struct)
%   method for Pq integration (1~4, 1 and 2 for 3D, 3 for 2D and 4 for 1D)
%
% When inputmode == 3(q for input), then Pqhkls and Pqhkls2 will be the
% integrated intensity for all particles. Otherwise Fqhkls or Pqhkls will
% be that of individual particles. So, SOF is included in the mode 3 only.
%
% particles.shape
%
% for cube
%   edge length of cube = 2*R, where R is particles.radius
%
%
% when anisotropic particles needs to be rotated....
% rotate pyramid 45degree in xy plane.
%x = qx;
%y = qy;
%qx = cos(xyangle)*x + sin(xyangle)*y;
%qy = -sin(xyangle)*x + cos(xyangle)*y;
% rotate pyramid 45degree in yz plane.
%y = qy;
%z = qz;
%qy = cos(yzangle)*y + sin(yzangle)*z;
%qz = -sin(yzangle)*y + cos(yzangle)*z;
%switch cellinfo.dim;
%    case 1
%        dim = 4;
%    case 2
%        dim = 3;
%    case 3
%        dim = 1;
%end

numparticle = numel(particles);
numP = zeros(numparticle, 1);
Np = 0;
isonlysphere = 1;
if ~iscell(particles)
    p = particles;
    particles = [];
    particles{1} = p;
end

for pindx=1:numparticle
    numP(pindx) = numel(particles{pindx}.position)/3*particles{pindx}.SOF;
    Np = Np + numP(pindx);   
    if ~strcmp(particles{pindx}.shape, 'sphere')
        isonlysphere = 0;
    end
end

%numPfrac = numP/Np;
numPfrac = numP;
mat = [];
if nargin < 4
    dim = 1;
end

if nargin>=3
    lvector = cellinfo.latticevectors;
    a1 = lvector(1, :);
    a2 = lvector(2, :);
    a3 = lvector(3, :);
    b1 = cross(a2, a3)/cellinfo.Vol;b1 = b1(:);
    b2 = cross(a3, a1)/cellinfo.Vol;b2 = b2(:);
    b3 = cross(a1, a2)/cellinfo.Vol;b3 = b3(:);
    B = [b1, b2, b3];
end

%typeofinput = 0;

if isstruct(hklsorq)
    hkls = hklsorq;
    typeofinput = 1;
    numq = numel(hkls);
else
    if size(hklsorq, 2) == 3
        hkls = [];
        qx = hklsorq(:,1);
        qy = hklsorq(:,2);
        qz = hklsorq(:,3);
        typeofinput = 2;
        numq = numel(qx);
        q = sqrt(qx.^2+qy.^2+qz.^2);
    else
        hkls = [];
        q = hklsorq;
        typeofinput = 3;
        numq = numel(q);
    end
end

switch typeofinput
    case 1 % when the input is the hkls.
        Fqhkls = zeros(numq, max([hkls.multiplicity]), numel(particles));
        Pqhkls = zeros(numq, numel(particles));
        Pqhkls2 = Pqhkls;
    case 2  % when the input is the q vector or qx,qy,qz
        if isonlysphere
            Fqhkls = zeros(numq, numel(particles{1}.radius), numel(particles));
        else
            Fqhkls = zeros(numq, 1, numel(particles));
        end
    case 3 % when the input is the total q.
        Fqhkls = [];
        Pqhkls = zeros(numq, numel(particles));
        Pqhkls2 = Pqhkls;
end
%oldshape = {};
for pshape = 1:numel(particles)
    particleshape = particles{pshape}.shape;
    isfromparticlemaker = 0;
    if isfield(particles{pshape}, 'edgelength')
        edL = particles{pshape}.edgelength;
        isfromparticlemaker = 1;
    end
    if isfromparticlemaker
        if isfield(particles{pshape}, 'edgelength_sig')
            sig = particles{pshape}.edgelength_sig;
        else
            sig = 0;
        end
    else
        if isfield(particles{pshape}, 'radius_sig')
            sig = particles{pshape}.radius_sig;
        else
            sig = 0;
        end
    end
    
    switch particleshape
        case 'sphere'
            R = particles{pshape}.radius;
            if isfield(particles{pshape}, 'radius_sig')
                sig = particles{pshape}.radius_sig;
            else
                sig = 0;
            end
            parameter = [R, sig];
            fffunctionname = 'sphere';
            if typeofinput ~= 3
                if ~isempty(hkls)
                    d = [hkls(:).D];
                    qn = 2*pi./d;
                    qn = qn(:);
                else
                    qn = q;
                end
                pt{1} = particles{pshape};
                [Fqtemp, Pqtemp] = SphFormFactorforSq(qn, pt);
                Fqhkls(1:size(Fqtemp,1), 1, pshape) = Fqtemp;
                Pqhkls(:, pshape) = Pqtemp;
                continue
            end
        case {'rhombicdodecahedron', 'rhombic dodecahedron'}
            if isfromparticlemaker
                R = edL/2;
                sig = sig/2;
            else
                R = particles{pshape}.radius;
            end
            parameter = R;
            fffunctionname = 'saxsrhombicdodecahedron';
            %Fqhkls(qindx, pp, pshape) = feval(fffunctionname, qx, qy, qz, parameter);

        case 'octahedron'
            if isfromparticlemaker
                R = edL/2;
                sig = sig/2;
            else
                R = particles{pshape}.radius;
            end
            parameter = R;
            % amplitude
            fffunctionname = 'saxsoctahedron';

        case 'truncatedoctahedron'
            tr = [];
            if isfromparticlemaker
                R = edL(1);
                t = edL(2);
                tr = t*R;
            else
                R = particles{pshape}.radius;
                if isfield(particles{pshape}, 'truncatedfraction')
                    t = particles{pshape}.truncatedfraction;
                    tr = t*R;
                end
            end
            parameter = [R, tr];

            %parameter = [R, f*R];
            % amplitude
            fffunctionname = 'saxstruncatedoctahedron';

        case 'truncatedcube'
            tr = [];
            if isfromparticlemaker
                R = edL(1);
                t = edL(2);
                tr = t*R;
            else
                R = particles{pshape}.radius;
                if isfield( particles{pshape}, 'truncatedfraction')
                    f = particles{pshape}.truncatedfraction;
                    tr = f*R;
                end
            end
            parameter = [R, tr];
            % amplitude
            fffunctionname = 'saxstruncatedcube';
            %Fqhkls(qindx, pp, pshape) = feval(fffunctionname, qx, qy, qz, parameter);
            % intensity
        case 'cube'
            if isfromparticlemaker
                R = edL(1);
                parameter = [R];
                if numel(edL) == 2;
                    H = edL(2);
                    parameter = [R, R, H];
                end
                parameter = parameter/2;
                sig = sig/2; % since radius becomes half.
            else
                R = particles{pshape}.radius(1);
                %H = particles{1}.height;
                parameter = R;
                if isfield(particles{pshape}, 'height')
                    H = particles{pshape}.height;
                    parameter = [R, R, H];
%                    obj = polyhedra('cuboid', R, H);
                end
            end
            % amplitude
            fffunctionname = 'saxscube';
            %Fqhkls(qindx, pp, pshape) = feval(fffunctionname, qx, qy, qz, parameter);
            % intensity
        case 'cylinder'
            H = [];
            if isfromparticlemaker
                R = edL(1)/2;
                if numel(edL) > 1
                    H = edL(2);
                end
                sig = sig/2;
            else
                R = particles{pshape}.radius(1);
                if isfield(particles{pshape}, 'height')
                    H = particles{pshape}.height;
                end
            end
            if isempty(H)
                error('field height is missing')
            end

%            R = particles{pshape}.radius;
%            H = particles{pshape}.height;
            parameter = [R, H];
            %[~, R] = rotate_around_vector([0,0,1], [1,0,0], 90);
            %t=[qx,qy,qz]*R;qx=t(1);qy=t(2);tz=t(3);
            fffunctionname = 'saxscylinder';
            %Fqhkls(qindx, pp, pshape) = feval(fffunctionname, qx, qy, qz, parameter);
        case 'concavecube'
            H = [];
            if isfromparticlemaker
                R = edL(1)/2;
                if numel(edL) > 1
                    H = edL(2);
                end
                sig = sig/2;
            else
                R = particles{pshape}.radius(1);
                if isfield(particles{pshape}, 'height')
                    H = particles{pshape}.height;
                end
            end
            if isempty(H)
                error('field height is missing')
            end

%            R = particles{pshape}.radius;
%            H = particles{pshape}.height;
            parameter = [edL(1), H];
            %[~, R] = rotate_around_vector([0,0,1], [1,0,0], 90);
            %t=[qx,qy,qz]*R;qx=t(1);qy=t(2);tz=t(3);
            fffunctionname = 'saxsconcavecube';
            %Fqhkls(qindx, pp, pshape) = feval(fffunctionname, qx, qy, qz, parameter);
        case 'convexcube'
            H = [];
            if isfromparticlemaker
                R = edL(1)/2;
                if numel(edL) > 1
                    H = edL(2);
                end
                sig = sig/2;
            else
                R = particles{pshape}.radius(1);
                if isfield(particles{pshape}, 'height')
                    H = particles{pshape}.height;
                end
            end
            if isempty(H)
                error('field height is missing')
            end

%            R = particles{pshape}.radius;
%            H = particles{pshape}.height;
            parameter = [edL(1), H];
            %[~, R] = rotate_around_vector([0,0,1], [1,0,0], 90);
            %t=[qx,qy,qz]*R;qx=t(1);qy=t(2);tz=t(3);
            fffunctionname = 'saxsconvexcube';
            %Fqhkls(qindx, pp, pshape) = feval(fffunctionname, qx, qy, qz, parameter);

        case 'THH'
            H = [];
            if isfromparticlemaker
                R = edL(1)/2;
                if numel(edL) > 1
                    H = edL(2)/2;
                end
                sig = sig/2;
            else
                if isfield(particles{pshape}, 'radius')
                    R = particles{pshape}.radius(1);
                end
                if isfield(particles{pshape}, 'height')
                    H = particles{pshape}.height;
                end
            end                    
            if isempty(H)
                H = R;
            end

            parameter = [R, H];
            %[~, R] = rotate_around_vector([0,0,1], [1,0,0], 90);
            %t=[qx,qy,qz]*R;qx=t(1);qy=t(2);tz=t(3);
            fffunctionname = 'saxsTHH';
            %Fqhkls(qindx, pp, pshape) = feval(fffunctionname, qx, qy, qz, parameter);
    end
% size distribution effect.....
    if sig~= 0
        [nr, xr] = schultzdist99(parameter(1), sig, 10);nr = nr/sum(nr);
    else
        xr = parameter(1); nr = 1;
    end
% =========================================

    if typeofinput ~= 3
        switch typeofinput
            case 1 % When hkls is input
                for qindx = 1:numel(hkls)
                    %for pp=1:numel(hkls(qindx).HKLs(1,:))
                    HKLv = hkls(qindx).HKLs;
                    zeroindx = find((HKLv(1,:)==0) & (HKLv(2,:)==0) & (HKLv(3,:)==0));
                    if ~isempty(zeroindx)
                        HKLv(:,zeroindx) = [];
                    end
                    %h = hkls(qindx).HKLs(1,:); % row vector
                    %k = hkls(qindx).HKLs(2,:);
                    %l = hkls(qindx).HKLs(3,:);
                    h = HKLv(1,:);
                    k = HKLv(2,:);
                    l = HKLv(3,:);
                %h = h(:);k=k(:);l=l(:);
                    %qv = 2*pi*(h*b1+k*b2+l*b3);
                    qv = 2*pi*B*[h;k;l];
                    qx = qv(1, :);
                    qy = qv(2, :);
                    qz = qv(3, :);

                    if ~isempty(mat)
                        temp= mat*[qx;qy;qz];
                        qx = temp(1, :);
                        qy = temp(2, :);
                        qz = temp(3, :);
                    end
                    qx = qx(:);
                    qy = qy(:);
                    qz = qz(:);

                    q_transform
                    Fq = calFq(fffunctionname, qxn, qyn, qzn, parameter, xr, nr);
                    Fq = Fq*particles{pshape}.rho;
                    %for kkk=1:numel(Fq);
                    %    Fqhkls(qindx, kkk, pshape) = Fq(kkk);
                    %end
                    if numel(particles)==1
                        Fqhkls(qindx, 1:numel(Fq)) = Fq';
                    else
                        Fqhkls(qindx, 1:size(Fq,1), pshape) = Fq';
                    end
                    Q = sqrt(qx.^2+qy.^2+qz.^2);
                    if isfield(particles{pshape}, 'Rotmat')
                        Rm{1} = particles{pshape}.Rotmat;
                    else
                        Rm{1} = eye(3);
                    end
                    [Pq, Pq2] = calPq(fffunctionname, Q(1), parameter, particles{pshape}, xr, nr, dim, Rm);
                    Pq = Pq*particles{pshape}.rho;
                    %Pq2 = Pq2*particles{pshape}.rho;
                    Pqhkls(qindx, pshape) = Pq;
                    %Pqhkls2(qindx, pshape) = Pq2;
                    %end
                end

            case 2
                q_transform
                Fq = calFq(fffunctionname, qxn, qyn, qzn, parameter, xr, nr);
                Fq = Fq*particles{pshape}.rho;
                if numel(particles) ==1
                    Fqhkls = Fq;
                else
                    Fqhkls(:, pshape) = Fq;
                end
                %Pq = calPq(fffunctionname, q, parameter, particles{pshape}, xr, nr);
                %Pqhkls(:, pshape) = Pq;
        
        end

    else
        fname{pshape} = fffunctionname;
        if isfield(particles{pshape}, 'Rotmat')
            Rm{pshape} = particles{pshape}.Rotmat;
        else
            Rm{pshape} = eye(3);
        end
        
        if ~strcmp(fffunctionname, 'sphere')
            n{pshape} = nr*numPfrac(pshape)*particles{pshape}.rho;
            x{pshape} = xr(:)/parameter(1)*parameter;
        else
            n{pshape} = numPfrac(pshape)*particles{pshape}.rho;
            x{pshape} = parameter;
        end
    end

end
if typeofinput == 3
    [Pq, Pq2] = calPqMultiparticle(fname, q, x, n, dim, Rm);
    Pqhkls = Pq;
    Pqhkls2 = Pq2;
end


try
    assignin('base', 'Fqhkls', Fqhkls)
    assignin('base', 'Pqhkls', Pqhkls)
end
    function q_transform
            if isfield(particles{pshape}, 'Rotmat')
                Rmat = particles{pshape}.Rotmat;
                %Q' = inv(mat) * [qx; qy; qz];
                Q = [qx(:), qy(:), qz(:)] * Rmat;% for rotation matrix R^T = R^-1.
                qxn = Q(:,1);
                qyn = Q(:,2);
                qzn = Q(:,3);
            else
                qxn = qx;
                qyn = qy;
                qzn = qz;
            end

            q = sqrt(qxn.^2+qyn.^2+qzn.^2);
            % Round q vector -- some function such as saxsTHH produce wrong
            % value when q is not rounded especially complex part of F.
            %qx = round(qx*1000000)/1000000;
            %qy = round(qy*1000000)/1000000;
            %qz = round(qz*1000000)/1000000;
    end
end
    
function Fq = calFq(fffunctionname, qx, qy, qz, parameter, xr, nr)
    if numel(xr) == 1
        Fq = feval(fffunctionname, qx, qy, qz, parameter);
    else
        %try % if the form factor can take multiple number of particles.
        %    p = xr(:)*parameter/parameter(1);
        %    tempint = feval(fffunctionname, qx, qy, qz, p);
        %    Fq = sum(repmat(nr, numel(qx), 1).*tempint, 2);
        %catch % Otherwise calculate one by one and add.
            Fq = zeros(numel(qx), 1);
            for isize=1:numel(nr)
                tempint = feval(fffunctionname, qx, qy, qz, parameter*xr(isize)/parameter(1));
                tempint(isnan(tempint)) = 0;
                Fq = Fq + tempint*nr(isize);
            end
        %end
        
    end
end

function [Pq, Pq2] = calPq(fffunctionname, q, parameter, obj, xr, nr, dim, Rm)
    if ~isfield(obj, 'Pq')
        if numel(xr) == 1
            nr = 1;
        end
        if nargout == 1
            Pq = saxs_average(q, fffunctionname, xr(:)/parameter(1)*parameter, nr, dim, Rm);
        else
            [Pq, Pq2] = saxs_average(q, fffunctionname, xr(:)/parameter(1)*parameter, nr, dim, Rm);
        end
%        else
%
%            Pq = zeros(numel(q),1);
%            for isize=1:numel(nr)
%                tempint = saxs_average(q, fffunctionname, parameter*xr(isize)/parameter(1));
%                Pq = Pq + tempint*nr(isize);
%            end
%        end
    else
        P = obj.Pq;
        Pq = interp1(P(:,1),P(:,2),q);
    end
end

function [Pq, Pq2] = calPqMultiparticle(fffunctionname, q, xr, nr, dim, Rm)

        [Pq, Pq2] = saxs_average(q, fffunctionname, xr, nr, dim, Rm);
%        else
%
%            Pq = zeros(numel(q),1);
%            for isize=1:numel(nr)
%                tempint = saxs_average(q, fffunctionname, parameter*xr(isize)/parameter(1));
%                Pq = Pq + tempint*nr(isize);
%            end
%        end
end
