function angle = q2angle(q, lambda, qvector, ai)
% function for converting q to angle(2th) 
% when q is a scalar.
%   angle = q2angle(q, lambda)
%       angle = rad2deg(asin(q*lambda/4/pi)*2)
% when q is [q, azimangle] % for powder diffraction
%   angle = q2angle(q, lambda, 1)
%       This covert [q, azimangle] into [tthf, azimangle]
% when q is [q, azimangle] % for powder diffraction
%   angle = q2angle(q, lambda, -1) % for powder diffraction
%       This covert [q, azimangle] into [tthf, af]
% when q is [qxy, qz] % for GI diffraction
%   angle = q2angle(q, lambda, 2)
%       This covert [qxy, qz] into [tthf, af]
% when q is [qx, qy, qz] % for GI diffraction
%   angle = q2angle(q, lambda, 3)
%       This covert [qx, qy, qz] into [tthf, af]
% if the geometry is a grazing incidence.
% use like this, angle = q2angle(q, lambda, 2 or 3, ai)      
%
% if you use polar coordinate.
% use like this, angle = q2angle([q, azimuthal angle], lambda, 1)      
% 
% default lambda is 1.5406
if nargin < 4
    ai = 0;
end
if nargin < 3
   qvector = 0;
end
if nargin < 2 
    lambda = 1.5406;
end

if qvector == 0;
    angle = rad2deg(asin(q*lambda/4/pi)*2);
elseif qvector == 1;   % convert [q, azimangle] to [tth, azimangle];
    qv = q(:,1);
    azimangle = q(:,2);
    ang = q2angle(qv, lambda);
    angle = [ang, azimangle];
elseif qvector == -1;   % convert [q, azimangle] to cartesian coordinate angles;
    error('Powder diffraction should not be described by cartesian angle coordinates')
    qv = q(:,1);
    azimangle = q(:,2);
    angle = 2*asin(lambda*[qv.*cos(azimangle), qv.*sin(azimangle)]/(4*pi))*180/pi;
    %[ang1, ang2] = pol2cart(azimangle, ang);
    %angle = [ang1(:), ang2(:)];
    %angle = [ang, azimangle];
elseif qvector == 2;   
    % convert [qxy, qz] to cartesian coordinate angles;
    ai = deg2rad(ai);
    qxy = q(:,1);
    qz = q(:,2);
    af = asin(qz*lambda/2/pi-sin(ai)); % 2af. 
    tthf = sign(qxy).*real(acos((cos(ai).^2 + cos(af).^2-(qxy*lambda/(2*pi)).^2)./(2*cos(ai)*cos(af))));
    ind = cellfun(@isreal, num2cell(tth));tthf(~ind) = NaN;
    angle = [rad2deg(tthf), rad2deg(af)];
    %angle = rad2deg(angle);
elseif qvector == 3;  % convert [qx, qy, qz] to cartesian coordinate angles;
                      % in this case, Ewald sphere is flat onto CCD.
                      % thus qx component is just neglected. 
    ai = deg2rad(ai);                      
    qy = q(:,2);
    qz = q(:,3);
    if qz*lambda/2/pi > pi/2
        disp('no')
    end
    af = asin(qz*lambda/2/pi-sin(ai)); % 2af. 
    if qz*lambda/2/pi-sin(ai) > pi/2
        disp('no')
    end
    tthf = asin(qy*lambda/(2*pi)./cos(af));
    
    angle = [rad2deg(tthf), rad2deg(af)];
    %angle = rad2deg(angle);
end