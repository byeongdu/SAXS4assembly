function d = dspacing(hkl, cellinfo, latticetype)
% d = dspacing(hkl, cellinfo)
% d = dspacing(hkl, cellinfo, latticetype)
% latticetype (string): cubic, tetragonal, hexagonal, orthorhombic,
% monoclinic, triclinic, trigonal(rhombohedral)--> hexagonal indexing
% from Culity book.
% hkl is an array of hkls for instance [1, 1, 1; 1, 2, 1; ...; 10, 1, 10];
% cellinfo is a structure containing, A, B, C, alpha, beta, gamma as
% fields.
%
% 2/11/2013
% Byeondu Lee

if isstruct(cellinfo)
    a = cellinfo.A;
    b = cellinfo.B;
    c = cellinfo.C;

    al = cellinfo.alpha*pi/180;
    be = cellinfo.beta*pi/180;
    ga = cellinfo.gamma*pi/180;
else
    a = cellinfo(1);
    b = cellinfo(2);
    c = cellinfo(3);
    al = cellinfo(4)*pi/180;
    be = cellinfo(5)*pi/180;
    ga = cellinfo(6)*pi/180;
end

h = hkl(:,1);
k = hkl(:,2);
l = hkl(:,3);

if a~=0 & b~=0 & c~= 0
    if isfield(cellinfo, 'Atensor')
        ds = sqrt(cellinfo.Atensor(1)*h.^2+...
            cellinfo.Atensor(2)*k.^2+...
            cellinfo.Atensor(3)*l.^2+...
            cellinfo.Atensor(4)*h.*k+...
            cellinfo.Atensor(5)*h.*l+...
            cellinfo.Atensor(6)*k.*l);
        d = 1./ds;
        return
    end
end
if nargin<3
    latticetype = 'triclinic';
end

cellp = [a~=0, b~=0, c~=0];
switch sum(cellp)
    case 1
        celldim = 1;
    case 2
        celldim = 2;
    case 3
        celldim = 3;
end

% if celldim == 1
%     latticepara = [a, b, c];
%     a = latticepara(cellp);
%     d = a./sqrt(h.^2+k.^2+l.^2);
%     return
% end


if a == 0
    a = 1E-9;
end
if b == 0
    b = 1E-9;
end
if c == 0
    c = 1E-9;
end

try
    if strcmp(strtrim(latticetype), 'trigonal')
        latticetype = 'hexagonal';% in my software, trigonal is indexed as hexagonal.
    end
    if strcmp(strtrim(latticetype), 'trigonal(R)')
        latticetype = 'trigonal2';% in my software, trigonal is indexed as hexagonal.
    end
catch
    latticetype = 'triclinic';
end
switch strtrim(latticetype)
    case 'cubic'
        d = a./sqrt(h.^2+k.^2+l.^2);
    case 'tetragonal'
        d = 1./sqrt((h.^2+k.^2)/a^2+l.^2/c^2);
    case 'hexagonal' 
%        if c~=0
%            d = 1./sqrt(4/3*(h.^2+h.*k+k.^2)/a^2+l.^2/c^2);
%        else
            d = 1./sqrt(4/3*(h.^2+h.*k+k.^2)/a^2);
%        end
    case 'orthorhombic'
        d = 1./sqrt(h.^2/a^2+k.^2/b^2+l.^2/c^2);
    case 'trigonal2'
       d = a*sqrt(1-3*cos(al)^2+2*cos(al)^3)./sqrt((h.^2+k.^2+l.^2)*sin(al)^2+...
           2*(h.*k+k.*l+h.*l)*(cos(al)^2-cos(al)));
    case 'monoclinic'
        d = sin(be)./sqrt((h/a).^2+(k*sin(be)/b).^2+(l/c).^2-2*h.*l*cos(be)/(a*c));
    case 'triclinic'
        V = a*b*c*sqrt(1-cos(al)^2-cos(be)^2-cos(ga)^2+2*cos(al)*cos(be)*cos(ga));
        S11 = (b*c*sin(al))^2;
        S22 = (a*c*sin(be))^2;
        S33 = (a*b*sin(ga))^2;
        S12 = a*b*c^2*(cos(al)*cos(be)-cos(ga));
        S23 = a^2*b*c*(cos(be)*cos(ga)-cos(al));
        S13 = a*b^2*c*(cos(ga)*cos(al)-cos(be));
        d = V./sqrt(S11*h.^2+S22*k.^2+S33*l.^2+...
            2*(S12*h.*k+S23*k.*l+S13*h.*l));
end
d(d<1E-5) = 0;