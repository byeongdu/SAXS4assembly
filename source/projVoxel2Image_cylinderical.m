function img = projVoxel2Image_cylinderical(Vox)
% in
Vox(isnan(Vox)) = 0;
sz = size(Vox);
imgsize = [round(sqrt(sum(sz(1:2).^2))), sz(3)];
[X, Y, ~] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
cenx = (sz(1)+1)/2;
ceny = (sz(2)+1)/2;
X = X-cenx;
Y = Y-ceny;
%Z = Z-sz(3)/2;
xy = sqrt(X.^2+Y.^2);%+sqrt((sz(1)/2)^2+(sz(2)/2)^2);
if imgsize(1)/2 == fix(imgsize(1)/2)
    xy_array = linspace(0, fix(max(xy(:,:,1), [], 'all'))+1, imgsize(1)/2+1);
    arraytype = 'even';
else
    xy_array = linspace(-0.5, fix(max(xy(:,:,1), [], 'all'))+1, (imgsize(1)-1)/2+2);
    arraytype = 'odd';
end
%xy_int = fix(xy);
%xy_frac = xy-xy_int;
%t = xy_int>0;
%xy_int(~t) = 170;
%ind = sub2ind(imgsize, xy_int, Z);
img = zeros(imgsize);%img = img';
%img2 = zeros(imgsize);
for i=1:sz(3)
    vd = Vox(:, :, i);
    [~, dd] = azimavg(vd, xy(:,:,1), xy_array, ones(size(xy(:,:,1))));
    if strcmp(arraytype, 'even')
        arr = [fliplr(dd'), dd'];
    else
        arr = [fliplr(dd'), dd(2:end)'];
    end
    img(:,i) = arr';
    % f = xy_frac(:, :, i);
    % %xy_int(xy_int==0) = 170;
    % [~, ~, c] = unique(xy_int(:, :, i));
    % img(:, i) = accumarray(c,vd(:).*(1-f(:)));
    % img(170, i) = 0;
    % xy_int(xy_int==170) = 0;
    % [~, ~, c] = unique(xy_int(:, :, i)+1);
    % img2 = accumarray(c,vd(:).*f(:));
    % img2(1, :) = 0;
%     [~, ~, c] = unique([xy_int(:), Z(:)], 'row', 'stable');
%     img = accumarray(c,Vox(:).*(1-xy_frac(:)));
%     img(170, :) = 0;
%     xy_int(xy_int==170) = 0;
%     [~, ~, c] = unique([xy_int(:)+1, Z(:)], 'row', 'stable');
%     img2 = accumarray(c,Vox(:).*xy_frac(:));
%     img2(1, :) = 0;
end
% img = img+img2;
% img = img';
%img = reshape(img, imgsize);