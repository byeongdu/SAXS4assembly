function test = guinier(varia, q)
% this for Guinier Function : guinier(variable, q)
% variable is composed of 3 parameter(eta, V, Rg)
% this function is valid in the region qRg < 1
%

%N = varia(1);
%eta = varia(2);
%V = varia(3);
%Rg = varia(4);

N = 1;
eta = varia(1);
V = varia(2);
Rg = varia(3);

test = N*eta^2*V^2.*exp(-Rg^2/3.*q.^2);