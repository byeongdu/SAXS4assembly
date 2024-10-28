function [qx, qy, qz] = angle2vq2(ai, af, thi, thf, lambda)
% function for converting angle(ai, af, theta) to vector q(qx, qy, qz)
% Grazing angle diffraction
% default lambda is 1.5406
if nargin < 3 
    lambda = 1.5406;
end
ai = deg2rad(ai);
af = deg2rad(af);
thetai = deg2rad(thi);
thetaf = deg2rad(thf);
tempqx = cos(af).*cos(thetaf) - cos(ai).*cos(thetai);
tempqy = cos(af).*sin(thetaf) + cos(ai).*sin(thetai);
tempqz = sin(af)+sin(ai);
qx = 2*pi/lambda*tempqx;qy = 2*pi/lambda*tempqy;qz = 2*pi/lambda*tempqz;
