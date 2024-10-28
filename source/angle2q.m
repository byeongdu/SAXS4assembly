function q = angle2q(angle, lambda)
% function for converting q to angle(2theta) : unit degree
% default lambda is 1.5406
if nargin < 2 
    lambda = 1.5406;
end
q = 4*pi/lambda*sin(deg2rad(angle/2));