function [Iq, Z, Sprime_q, hkls, Pprime_q] = strfactor_random(q, hkls, particles, cellinfo, NN, siga, factorX, peakw)
j = sqrt(-1);
hkls(end) = [];
isSphericalParticle = 1;

cellpara = [cellinfo.A, cellinfo.B, cellinfo.C];
dim_cell = numel(find(cellpara~=0));
dim_orientation = 3;

%if dim_cell == 2, dim_orientation can be either 2 or 1;
%if dim_cell == 1, dim_orientation should be 1;
if dim_orientation > dim_cell
    disp('Dimensionality of random orientation should be equal or less than the dimensionality of cell')
    dim_orientation = dim_cell;
    fprintf('The mfile strfactor.m sets the dimensionality of random orientation to %dD powder\n', dim_cell);
end
% voigt
% p = [ Amplitude Centre Gauss_width Lorz_width Background ]
% atomic scattering factor
if nargin < 6
    siga = 0.05;
    factorX = 10;
    peakw = 1/NN*0.2;
end

if nargin == 5
    gw = 0.01;
    lw = 0.02;
%    gw = (2*pi/NN^2*0.01);
%    lw = 2*gw;
else
    gw = 0.01;
    lw = 0.02;
end

% peak width
%lw = 1/NN*0.2;
lw = peakw;
gw = lw*0.5;


qmax = max(q);
numparticle = numel(particles);

nP = 0;
d = [hkls(:).D];
qdata = 2*pi./d;
qdata = qdata(:);
Fsp = zeros(numel(qdata), numparticle);

Np = 0;
Rad = zeros(numparticle, 1);
numP = zeros(numparticle, 1);
for pindx=1:numparticle
    numP(pindx) = numel(particles{pindx}.position)/3;
    Rad(pindx) = particles{pindx}.radius;
   Np = Np + numel(particles{pindx}.position)/3;
end
lvector = cellinfo.latticevectors;
a1 = lvector(1, :);
a2 = lvector(2, :);
a3 = lvector(3, :);
b1 = cross(a2, a3)/cellinfo.Vol;
b2 = cross(a3, a1)/cellinfo.Vol;
b3 = cross(a1, a2)/cellinfo.Vol;

NeedtoCalFq = 1;
if isfield(particles{pindx}, 'Fq') 
    NeedtoCalFq = 0;
end
NeedtoCalFq = 1;

[Fsp, Fsphkls] = SphFormFactorforSq(q, particles, hkls);

for pindx=1:numparticle
   Rd = particles{pindx}.radius;
   TfactorB = particles{pindx}.TfactorB;
   nP = numel(particles{pindx}.position)/3;
   %disp(sprintf('The number of particle with radius %0.3f is %d', Rd, nP))
   %SOF = particles{pindx}.SOF;
   if ~NeedtoCalFq
   %    Fsp(:,pindx) = interp1(q, particles{pindx}.Fq(:), qdata);
       Tfactor(:,pindx) = exp(-TfactorB*qdata.^2);
   else
   %    sp = sphereamp(qdata,Rd);
   %    Fsp(:,pindx) = sp(:);
       Tfactor(:,pindx) = exp(-TfactorB*qdata.^2);
   end
   %Fsp(:,pindx) = Fsp(:,pindx)*SOF;
end

if factorX <= 0
    factorX = nP;
end

Fhkl2 = zeros(size(d));
mat = particles{1}.mat';
am = [];
ppos = [];
nV2 = 0;
nV2SOF = 0;
for pindx=1:numparticle
    nP = numel(particles{pindx}.position)/3;
    ppos = [ppos; particles{pindx}.position*mat];
    am = [am; ones(nP,1)*pindx];
    SOF = particles{pindx}.SOF;
    nV2SOF = nV2SOF + nP*(SOF*4*pi/3*particles{pindx}.radius.^3).^2;
    nV2 = nV2 + nP*(4*pi/3*particles{pindx}.radius.^3).^2;
end


fhkl_random
assignin('base', 'hkls', hkls)

Z = zeros(size(q));

switch dim_orientation
    case 1
        omega = 1;
    case 2
        omega = 2*pi;
    case 3
        omega = 4*pi;
    otherwise
        1;
end

% Exact definition of beta(q) = |Fu(q)|^2/P'(q)
Fm = [];

%Pprime_q = SchultzMultiSphereFun(q, Rad, 0.0001, numP);
Pprime_q = zeros(size(q));

for pindx=1:numparticle
    SOF = particles{pindx}.SOF;
    Fm = [Fm, SOF*Fsp(:, pindx);];
%    Fm = [Fm, SOF*sphereamp(q,Rad(pindx));];
%    Fm = [Fm, sphereamp(q,Rad(pindx));];
    Pprime_q = Pprime_q + SOF*numP(pindx)*abs(Fsp(:, pindx)).^2;
%    Pprime_q = Pprime_q + SOF*numP(pindx)*spherePq(q, Rad(pindx));
%    Pprime_q = Pprime_q + numP(pindx)*spherePq(q, Rad(pindx));
end
Fu2_q = debyem(q, ppos(:,1), ppos(:,2), ppos(:,3), am, Fm);
betaq = Fu2_q./Pprime_q;

factor = 1/omega*(2*pi)^dim_cell/cellinfo.Vol;
Lorentzfactor = qdata.^(dim_orientation-1);
for i=1:numel(d)
    qd = Lorentzfactor(i);
    %peak = lw/(2*pi)./((q-qdata(i)).^2 + (lw/2)^2);
    peak = lorza(q, [factorX, qdata(i), lw, 0]);
    %peak = voigt(q, [1, qdata(i), gw, lw, 0])/(gw + gw*0.06);
    Z0 = Fhkl2rP(i)*peak/qd*factor;
    Z = Z + Z0;
    %Z2 = Z2 + Fhkl3rP(i)*peak/qd*factor/12;
end
%Z = Z - Pq(:)./Fq(:);
%Z = zeros(size(Z));
%Z = Z+1;
%Z = Z./Fu2_q;

if nargin > 4
    Gq = exp(-siga^2*(NN)^2*q.^2);
else
    Gq = 1;
end

% polydispersity of particles: beta(q) = exp(-sig^2*R^2*q^2)
%Rd = max(radius);
%sigR = 0.02;
%sigR = 0.08;
%beta = exp(-sigR^2*Rd^2*q.^2);
beta = 1; % because it is not polydisperse



Iq = 1+(Z - 1).*Gq;

%Sprime_q = betaq.*(Z-1).*Gq + 1;
Sprime_q = betaq;%.*(1+(Z-1).*Gq);

if nargin < 5
    Iq = Z;return
else
    if isempty(NN)
        return
    end
    if NN <= 0
        return
    end
end



function fhkl_random
nP = 0;
numRandom = 1000;
tnp = 0; %total number of particles
for i=1:numparticle
    tnp = tnp + size(particles{i}.position, 1);
end

%fhkl = zeros(tnp, 1);
%fhkla = zeros(numel(d), 1);
atomPOS = zeros(tnp, 3);
randv = zeros(numRandom, tnp);
fhkl = zeros(numRandom, tnp);

for i=1:numparticle
        SOF = particles{i}.SOF;
        atm = particles{i}.position;
        numatom = size(atm, 1);
        atm = particles{i}.position;
        for k=1:numatom
            nP = nP+1;
            if SOF>0
                randv(:, nP) = fix(rand(numRandom,1)+SOF);
            else
                randv(:, nP) = zeros(numRandom, 1);
            end
            atomPOS(nP, :) = atm(k, :);
        end
end

for qindx = 1:numel(d)
    if (qdata(qindx) > qmax)
        break
    end
    qm = qdata(qindx);
    nP = 0;
    for i=1:numparticle
        fhkl0 = Fsphkls(qindx,i);
%        fhkl0 = fhkl0*Tfactor(qindx, i);%temperature factor
        atm = particles{i}.position;
        numatom = size(atm, 1);
        for k=1:numatom
            nP = nP+1;
            fhkl(:,nP) = randv(:, nP)*fhkl0;
        end
    end
%    fhkl = fhkl';  %matrix for the atomic scattering factor; [100, nP] size.
    hkl = hkls(qindx).HKL;
    hkl = hkl(:);
    fhkla = fhkl*exp(2*pi*j*atomPOS*hkl);
    Ihkla = abs(fhkla).^2;
    Ihkla = sum(Ihkla)/numRandom;
    hkls(qindx).I = Ihkla;
    Fhkl2(qindx) = hkls(qindx).multiplicity.*Ihkla;

    for pindx=1:numparticle
        SOF = particles{pindx}.SOF;
        FF(pindx) = Fsphkls(qindx,pindx)*SOF;
    end

    Pq = debyem(qm, ppos(:,1), ppos(:,2), ppos(:,3), am, FF);
    %Pq2 = SchultzMultiSphereFun(qm, Rad, 0.1, numP);
    Fhkl2rP(qindx) = Fhkl2(qindx)/Pq; % the form factor is removed.    
end

end
end

