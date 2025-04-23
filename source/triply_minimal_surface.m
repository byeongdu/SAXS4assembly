function [f, p] = triply_minimal_surface(type, a, N)
if nargin < 2
    a = 1;
end
if nargin <3
    N = 2^6;
end
% Parameters you want to play with...
%a = 731.3433; % lattice constant.
%volf1 = 0.15;
%volf2 = 0.70;
%volf3 = 1 - volf1 - volf2;
a_G = 2*pi;
x = linspace(0, a, N); % periodiciy should be 2*pi.
y = linspace(0, a, N);
z = linspace(0, a, N);
% x = linspace(0, a, N+1); % periodiciy should be 2*pi.
% y = linspace(0, a, N+1);
% z = linspace(0, a, N+1);
% x(end) = [];
% y(end) = [];
% z(end) = [];
[X, Y, Z] = ndgrid(x, y, z);
switch type
    case 'P'
        % P-surface
        f = cos(X/a*a_G) + cos(Y/a*a_G) + cos(Z/a*a_G);
    case 'D'
        % D-surface
        f = cos(X/a*a_G).*cos(Y/a*a_G).*cos(Z/a*a_G) -...
            sin(X/a*a_G).*sin(Y/a*a_G).*sin(Z/a*a_G);
    case 'G'
        % G-surface
% x = linspace(-a/4, 5/4*a, 101); % periodiciy should be 2*pi.
% y = linspace(-a/4, 5/4*a, 101);
% z = linspace(-a/4, 5/4*a, 101);
% [X, Y, Z] = ndgrid(x, y, z);

        f = sin(X/a*a_G).*cos(Y/a*a_G) + sin(Y/a*a_G).*cos(Z/a*a_G) + ...
            sin(Z/a*a_G).*cos(X/a*a_G);
% p = isosurface(X,Y,Z,f,0);
% N_cm = length(p.faces);
% cm = zeros(N_cm, 3);
% for i=1:N_cm
%     cm(i, :) = mean(p.vertices(p.faces(i, :), :));
% end
    case 'I-WP'
         %I-WP surface
% a = 34; % lattice constant.
% a_G = 2*pi;
% N = 92;
% x = linspace(0, a, N); % periodiciy should be 2*pi.
% y = linspace(0, a, N);
% z = linspace(0, a, N);
% [X, Y, Z] = ndgrid(x, y, z);

        f = 2*(cos(X/a*a_G).*cos(Y/a*a_G) + cos(Z/a*a_G).*cos(X/a*a_G)+cos(Y/a*a_G).*cos(Z/a*a_G))...
            -(cos(2*X/a*a_G)+cos(2*Y/a*a_G)+cos(2*Z/a*a_G));
end
p = isosurface(X,Y,Z,f,0);
p.facecolor = 'r';
p.facealpha = 0.8;
p.edgecolor = 'none';
%figure
%patch('Faces', p.faces, 'Vertices', p.vertices, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
%camlight;
%lighting gouraud;
%camva(8);
% %% 3D unitcell
% thick1 = 0.022*volf1/0.1331*a; % matrix
% thick2 = 0.1739*volf2/0.7666*a; % shell 
% thick3 = 0.25*volf3/0.1003*a; % core
% Npix = 100;
% 
% N_cm = length(cm);
% x = linspace(0, a, Npix+1); % periodiciy should be 2*pi.
% x(end) = [];
% [X, Y, Z] = ndgrid(x, x, x);
% ph = gyroid_fillcell(cm, X(:), Y(:), Z(:), [thick1, thick2, thick3]);
% ph = reshape(ph, [Npix, Npix, Npix]);
% Nvox = Npix^3;
% disp('Done')
% for i=1:3
%     fprintf('Volume fraction of phase %i of the constructed unit cell is %0.3f\n', i, sum(ph(:)==i)/Nvox);
% end
% save GY.mat ph X Y Z a x y z