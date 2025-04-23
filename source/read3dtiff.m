function dd = read3dtiff(filename)
obj = Tiff(filename, 'r');
currentDirectory(obj);
dd = [];i = 1;
while 1
    dd(:,:,i) = obj.read;
    if lastDirectory(obj)
        break
    end
    nextDirectory(obj);
    i = i+1;
end
close(obj)