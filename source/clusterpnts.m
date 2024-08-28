function pnts = clusterpnts(pnts, ds)
% pnts can be XYZ or XY
% ds is the minimum distance.
if nargin<2
    ds = 2;
end
    k = 1;
    while k<size(pnts, 1)
        p0 = pnts(k, :);
        indx = cluster(pnts, p0, ds);
        try
        pnts(indx, :) = repmat(mean(pnts(indx, :),1), sum(indx), 1);
        catch
            pnts;
        end
        indx(k) = 0;
        pnts(indx, :) = [];
        k = k+1;
    end


    function ind = cluster(pnts, pnt, ds)
        d = sqrt(sum((pnt-pnts).^2, 2));
        ind = d < ds;
    end

end