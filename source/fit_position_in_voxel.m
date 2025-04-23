function [pos, Iq0, sig] = fit_position_in_voxel(Iq, X, Y, Z, i, j, k)
% [pos, Iq, sig] = fit_position_in_voxel(Vox, X, Y, Z, i, j, k)
% input: 
%   Vox : voxel
%   X : x coordinate of the vox
%   Y : y coordinate of the vox
%   Z : z coordinate of the vox
%   i, j, k: initial positions of the center
% output:
%   pos : fitted center position
%   Iq : integrated intensity (gaussain function)
%   sig : std of the 3d gaussian function
%

% when a matrix is given, it returns the peak position as an index.
    options = optimset('fminsearch');
    options = optimset(options, 'TolX',0.1E-6);
    options = optimset(options, 'MaxIter',5000);
    options = optimset(options, 'MaxFunEvals', 5000);
    dX = abs(X(2,1,1)-X(1,1,1));
    dY = abs(Y(1,2,1)-Y(1,1,1));
    if dX == 0
        dY = abs(Y(2,1,1)-Y(1,1,1));
        dX = abs(X(1,2,1)-X(1,1,1));
    end
    dZ = abs(Z(1,1,2)-Z(1,1,1));
    dv = dX*dY*dZ;
    mv = abs(sum(Iq(:))*dv);
    xl = [min(X(:)), max(X(:))]; % x range
    yl = [min(Y(:)), max(Y(:))];  % y range
    zl = [min(Z(:)), max(Z(:))];  % zrange
    %mindim = min([diff(xl), diff(yl), diff(zl)]);
    Iq = Iq(:);
    NLPstart = [mv, i, j, k, dX, dY, dZ, mv/100]; % I0, x0, y0, z0, sigma
    LB = [mv*0.1, xl(1), yl(1), zl(1), dX/10, dY/10, dZ/10, -mv*0.1];
    UB = [mv*10, xl(2), yl(2), zl(2), dX*10, dY*10, dZ*10, mv*0.1];
%     minpixeldistance = min([dX, dY, dZ]);
%     NLPstart = [mv, i, j, k, minpixeldistance/2]; % I0, x0, y0, z0, sigma
%     LB = [mv*0.1, xl(1), yl(1), zl(1), minpixeldistance/10];
%     UB = [mv*10, xl(2), yl(2), zl(2), minpixeldistance*2];
    pos = [X(:), Y(:), Z(:)];
    y = fminsearchcon(@(x) fitfunc(x, Iq, pos),NLPstart,LB,UB, [], [], [], options);
    pos = y(2:4);
    Iq0 = y(1);
    sig = y(5:7);
%%
function err = fitfunc(param, Iq, pos)
    % pos: mx3 matrix
    % Iq : mx1 matrix
    % param: 1x5 matrix
    cen = abs(param(2:4));
    sig = abs(param(5:end-1));
    back = param(end);
    ycal = param(1)*gaussian3d(pos(:,1), pos(:,2), pos(:,3), cen, sig)+back;
    err = chi_squared(Iq, ycal, 4);