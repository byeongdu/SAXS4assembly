function run_inflip(inflipfile)
[~, sg, fn, ~, data] = loadinflip(inflipfile);
hkl = data(:, 1:3);
Iq = data(:,4);
%[gk, chi2, err, gk0] = solve_densitymap(hkl,Iq, 'display',true);
[gk, chi2, err, gk0] = solve_densitymap(hkl,Iq, 'spacegroup', sg, 'display',true);

gk(end, :, :)=[];
gk(:, end, :)=[];
gk(:,:, end) =[];
WriteMRC(round(gk/max(gk(:))*255), 3, fn, 1);
drawccp4(inflipfile)