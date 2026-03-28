function y = pvza(x, p)
% pvza : pseudo-Voigt Area (1D)
% y = pvza(x, p)
% p = [Area, Centre, FWHM, eta, Background]
%
% Pseudo-Voigt: pv(x) = eta*L(x) + (1-eta)*G(x)
%   eta = 0 : pure Gaussian (shorter tails)
%   eta = 1 : pure Lorentzian (longer tails)
%
% Both components share the same FWHM (p(3)).
% The Lorentzian uses FWHM directly (as in lorza).
% The Gaussian converts FWHM -> sigma = FWHM / (2*sqrt(2*log(2))).
%
% Area is conserved: integral of y-Background = p(1) for any eta.
%
% Reference: Thompson, Cox & Hastings (1987); Dinnebier & Scardi (2021).

eta = p(4);
bg  = p(5);

% Lorentzian component (lorza uses p(3) as FWHM)
y_L = lorza(x, [p(1), p(2), p(3), 0]);

% Gaussian component: convert FWHM to sigma
sigma = p(3) / (2 * sqrt(2 * log(2)));
y_G = gaussa(x, [p(1), p(2), sigma, 0]);

y = bg + eta .* y_L + (1 - eta) .* y_G;
