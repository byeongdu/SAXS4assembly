function ct = celcon4rcplatvec(as, bs, cs)
% function ct = celcon4rcplatvec(as, bs, cs)
    rcpV = dot(as, cross(bs, cs));

    a = (2*pi)*cross(bs, cs)/rcpV;
    b = (2*pi)*cross(cs, as)/rcpV;
    c = (2*pi)*cross(as, bs)/rcpV;
    ct.recilatticevectors = [as; bs; cs];
    ct.recimat = ct.recilatticevectors';
    ct.latticevectors = [a; b; c];
    ct.mat = [a',b',c'];
    A = norm(a);
    B = norm(b);
    C = norm(c);
    alpha = angle2vect2(b, c)*180/pi;
    beta = angle2vect2(a, c)*180/pi;
    gamma = angle2vect2(a, b)*180/pi;
    ct.A = A;
    ct.B = B;
    ct.C = C;
    ct.alpha = alpha;
    ct.beta = beta;
    ct.gamma = gamma;
    RCP = [A, B, C, alpha, beta, gamma];
    
    RAD = 180/pi;
    AL=cos(alpha/RAD);
    BE=cos(beta/RAD);
    GA=cos(gamma/RAD);

    VOL=A*B*C*sqrt(1.0-AL^2-BE^2-GA^2+2.0*AL*BE*GA);
    ct.Vol = VOL;    

    % calculate orientation matrix
    cp = [A,B,C,alpha,beta,gamma];
    t = celcon(cp);Bmat=t.recimat;
    UB = ct.recimat;
    U = UB*inv(Bmat);
    ct.U = U;

    
    % https://gsas-ii.readthedocs.io/en/latest/GSASIIutil.html
% The “A tensor” terms are defined as A=(G11,G22,G33,2G12,2G13,2G23) 
% and A can be used in this fashion: d∗=A0h2+A1k2+A2l2+A3hk+A4hl+A5kl
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−√, 
% where d is the d-spacing, and d∗ is the reciprocal lattice spacing, 
% Q=2πd∗=2π/d. 
% g is the lattice tensor.
g = [A^2, A*B*GA, A*C*BE; 
    A*B*GA, B^2, B*C*AL; 
    A*C*BE, B*C*AL, C^2];
G = inv(g); % G is the reciprocal lattice tensor.
% A tensor is defined as
Atensor = [G(1,1), G(2,2),G(3,3),2*G(1,2),2*G(1,3),2*G(2,3)];

ct.g = g;
ct.G = G;
ct.Atensor = Atensor;
ct.cell = RCP;    