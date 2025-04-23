function patchvoxels(coordinates, map, map_range)
if max(map(:))>255
    map = fix(map/max(map(:))*255);
end
t = find(map>map_range(1) & map<map_range(2));
cm = jet(diff(map_range));
figure;
colormap(jet)
sprintf('Total number of points to draw: %i\n', numel(t))
for i=1:numel(t)
    voxel(coordinates(t(i), :), cm(map(t(i))-map_range(1), :),0.1);
end
colorbar;
axis image;