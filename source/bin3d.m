function data = bin3d(data, res)
if nargin<2
    res = 2;
end
if ndims ==3
    x = 1:size(data, 1);
    x = ceil(x/res);
    y = 1:size(data, 2);
    y = ceil(y/res);
    z = 1:size(data, 3);
    z = ceil(z/res);
    
    [X, Y, Z] = ndgrid(x, y, z);
    [~,~,IC] = unique_m([X(:), Y(:), Z(:)]);
elseif ndims ==2
    x = 1:size(data, 1);
    x = ceil(x/res);
    y = 1:size(data, 2);
    y = ceil(y/res);
    
    [X, Y, Z] = ndgrid(x, y, z);
    [~,~,IC] = unique_m([X(:), Y(:), Z(:)]);
end
dat = accumarray(IC, data(:));
h = histcounts(IC, 1:(max(IC)+1));
dat = dat./h(:);
data = reshape(dat,[max(x), max(y), max(z)]);