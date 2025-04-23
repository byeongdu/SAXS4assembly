function ct = celcon4latvec(a, b, c)

is2d = false;
if nargin == 2
    is2d = true;
else
    if sum(c.^2)==0
        is2d = true;
    end
end
if is2d
    c = [0, 0, 1];
end
A = norm(a);
B = norm(b);
if is2d
    C = 0;
else
    C = norm(c);
end

if is2d
    alpha = 0;
    beta = 0;
else
    alpha = angle2vect2(b, c)*180/pi;
    beta = angle2vect2(a, c)*180/pi;
end
gamma = angle2vect2(a, b)*180/pi;
RCP = [A, B, C, alpha, beta, gamma];

V = dot(a, cross(b, c));

aR = 2*pi*cross(b, c)/V;
bR = 2*pi*cross(c, a)/V;
cR = 2*pi*cross(a, b)/V;
if is2d
    cR = [0, 0, 0];
    c = [0, 0, 0];
end
ct.latticevectors = [a; b; c];
ct.mat = ct.latticevectors';
ct.recilatticevectors = [aR; bR; cR];
ct.recimat = [aR',bR',cR'];
ct.A = A;
ct.B = B;
ct.C = C;
ct.alpha = alpha;
ct.beta = beta;
ct.gamma = gamma;
RAD = 180/pi;
AL=cos(alpha/RAD);
BE=cos(beta/RAD);
GA=cos(gamma/RAD);

%VOL=A*B*C*sqrt(1.0-AL^2-BE^2-GA^2+2.0*AL*BE*GA);
VOL = V;
ct.Vol = VOL;

% calculate orientation matrix
cp = [A,B,C,alpha,beta,gamma];
t = celcon(cp);Bmat=t.recimat;
UB = ct.recimat;
if is2d
    invB = inv(Bmat(1:2, 1:2));
    invB = padarray(invB, [1, 1], 0, 'post');
else
    invB = inv(Bmat);
end
U = UB*invB;
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
if is2d
    invg = inv(g(1:2, 1:2));
    G = padarray(invg, [1, 1], 0, 'post');
else
    G = inv(g); % G is the reciprocal lattice tensor.
end
% A tensor is defined as
Atensor = [G(1,1), G(2,2),G(3,3),2*G(1,2),2*G(1,3),2*G(2,3)];

ct.g = g;
ct.G = G;
ct.Atensor = Atensor;
ct.cell = RCP;