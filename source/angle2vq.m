function q = angle2vq(alpha, theta, lambda)
% function for converting angle(ai, af, theta) to vector q(qx, qy, qz)
% Grazing angle diffraction
% default lambda is 1.5406
if nargin < 3 
    lambda = 1.5406;
end
ai = deg2rad(alpha(:, 1));
af = deg2rad(alpha(:, 2));
thetai = deg2rad(theta(:, 1));
thetaf = deg2rad(theta(:, 2));
tempq = [cos(af).*cos(thetaf) - cos(ai).*cos(thetai), cos(af).*sin(thetaf) + cos(ai).*sin(thetai), sin(af)+sin(ai)];
q = 2*pi/lambda*tempq;