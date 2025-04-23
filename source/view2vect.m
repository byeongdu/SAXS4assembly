function [as, Rmat] = view2vect(fh)
if nargin < 1
    fh = gcf;
end
ax = findobj(fh, 'type', 'axes');
[az, el] = view(ax);
Rmat = viewmtx(az, el);
vz = sin(el*pi/180);
vx = cos(el*pi/180)*sin(az*pi/180);
vy = -cos(el*pi/180)*cos(az*pi/180);
as = [vx, vy, vz];
as = as/norm(as);