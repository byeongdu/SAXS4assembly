function I = pvza3D(x, y, z, p)
% pvza3D : pseudo-Voigt Area (3D)
% I = pvza3D(x, y, z, p)
% p = [Area, x0, y0, z0, Width, eta]
%
% Pseudo-Voigt: pv = eta*L + (1-eta)*G
%   eta = 0 : pure 3D Gaussian (Gaussian3d, shorter tails)
%   eta = 1 : pure 3D Lorentzian (lorza3D, longer tails)
%
% Width (p(5)) plays the role of:
%   Lorentzian: gamma (half-width parameter, as in lorza3D)
%   Gaussian:   sigma (standard deviation, same numeric value)
%
% Reference: Thompson, Cox & Hastings (1987); Dinnebier & Scardi (2021).

eta = p(6);

% Lorentzian component
I_L = lorza3D(x, y, z, p(1:5));

% Gaussian component (gaussian3d returns unit-area pdf; scale by Area)
sigma = p(5);
I_G = p(1) * gaussian3d(x, y, z, [p(2), p(3), p(4)], sigma);

I = eta .* I_L + (1 - eta) .* I_G;
