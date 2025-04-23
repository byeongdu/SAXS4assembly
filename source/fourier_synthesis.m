function rho = fourier_synthesis(hkl, Fhkl, varargin)
% rho = fourier_synthesis(hkl, Fhkl, N)
% rho = fourier_synthesis(hkl, Fhkl, N, cellinfo)
% rho = fourier_synthesis(hkl, Fhkl, fracX, fracY, fracZ)
% rho = fourier_synthesis(hkl, Fhkl, fracX, fracY, fracZ, cellinfo)
% rho : density of the unit cell (range from 0 to N-1/N)
% if cellinfo is not given, the volume of the cell will be 1.
%
% Byeongdu Lee
% Fhkl =Fhkl*exp(iϕ_hkl)
% d(ha*)d(kb*)d(lc*)=a*b*c*dhdkdl=Vcell^*∗dhdkdl= 1/Vcell*dhdkdl
% h,k,l∈ ⇒ dh=dk=dl=Δh=Δk=Δl=1
% ρ(x,y,z)= 1/Vcell Sigma_-inf^inf|h Sigma_-inf^inf|k Sigma_-inf^inf|l
% (Fhkl*cos(ϕhkl −2π(hx+ky+lz))
N = 1;
Vcell = 1;
frac_X = [];
frac_Y = [];
frac_Z = [];
if numel(varargin) == 1
    N = varargin{1};
end
if numel(varargin) == 2
    cellinfo = varargin{2};
    N = varargin{1};
    Vcell = cellinfo.Vol;
end
if numel(varargin) == 3
    frac_X = varargin{1};
    frac_Y = varargin{2};
    frac_Z = varargin{3};
    N = size(frac_X);
    sz = N;
end
if numel(varargin) == 4
    frac_X = varargin{1};
    frac_Y = varargin{2};
    frac_Z = varargin{3};
    cellinfo = varargin{4};
    Vcell = cellinfo.Vol;
    N = size(frac_X);
    sz = N;
end
ndN = numel(N);
% if frac_X is not provided, generate.
if isempty(frac_Z)
    if ndN == 1
        if N > 1
            x = linspace(0, 1, N+1);
            y = linspace(0, 1, N+1);
            z = linspace(0, 1, N+1);
            x(end) = [];y(end)=[];z(end)=[];
            [frac_X, frac_Y, frac_Z] = ndgrid(x, y, z);
            sz = size(frac_X);
    %         frac_X = frac_X(:)/sz(2);
    %         frac_Y = frac_Y(:)/sz(1);
    %         frac_Z = frac_Z(:)/sz(3);
        end
    end
    if ndN == 2
        x = linspace(0, 1, N(1)+1);
        y = linspace(0, 1, N(2)+1);
        x(end) = [];y(end)=[];
        z = 0;
        [frac_X, frac_Y, frac_Z] = ndgrid(x, y, z);
        sz = size(frac_X);
    end
end
%V = cellinfo.Vol;
switch ndN
    case 1
        rho = zeros(size(frac_X));
    case 2
        rho = zeros(N(1), N(2));
    case 3
        rho = zeros(N(1), N(2), N(3));
end
if size(hkl, 2) == 2
    hkl(:,3) = zeros(1, size(hkl, 1));
    frac_Z = 0;
end

% F(hkl) = sum_x sum_y sum_z rho(x,y,z)*exp(-j*2*pi*(hx+ky+lz))
% rho(x, y, z) = sum_h sum_k sum_l |F(hkl)|exp(j*2*pi*(hx+ky+lz))
% Thus, phase angle of F(hkl) is opposite sign of the phase angle in
% conventional crystallographic definition.
% If, the phaseangle is defined as angle(F(hkl)), then
% rho(x, y, z) = sum_h sum_k sum_l |F(hkl)|exp(-j*phaseangle)
%
% Therefore,
% rho(x, y, z) = sum_h sum_k sum_l |F(hkl)|*exp(2*pi*(hx+ky+lz))
% rho(x, y, z) = sum_h sum_k sum_l |F(hkl)|*cos(2*pi*(hx+ky+lz) + phaseangle)

for i=1:size(hkl, 1)
    h = hkl(i, 1);
    k = hkl(i, 2);
    l = hkl(i, 3);
    if size(Fhkl, 2) == 1
        % frac = [0,0,0] --> exp(2*pi*j*(h*0+k*0+l*0)) == 1
        rhot = Fhkl(i)*exp(2*pi*sqrt(-1)*(h*frac_X+k*frac_Y+l*frac_Z));
    else
        rhot = Fhkl(i, 1)*cos(2*pi*(h*frac_X+k*frac_Y+l*frac_Z)+Fhkl(i, 2));
    end
    rho = rho + rhot;
end
rho = reshape(rho, sz)/Vcell;