function Fq = atomicscatteringfactor(q, atomN, charge)
try
    load atomicscatteringfactor.mat
catch
    % when failed to load the mat file, the information will be loaded from
    % the table file.
    fid = fopen('atomicscatteringfactor_coefficient.txt', 'r');
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    C = textscan(fid, '%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f');
    fclose(fid);
    atmN = zeros(size(C{1}));
    a = atom;
    for i=1:numel(C{1})
        t = findcellstr(a.Sym_, C{1}(i)); 
        atmN(i)=t;
    end
    [~, neutralatom, allatoms] = unique(atmN);
    t = cellfun(@(x) numel(x), C{2})-1;
    k = find(t==1);
    m = cell2mat(C{2}(k));
    m = [m(:,2), m(:,1)];
    for i=1:numel(k)
        C{2}(k(i)) = {m(i, :)};
    end
    charge_atoms = cellfun(@(x) str2double(x), C{2});
    charge_atoms(isnan(charge_atoms)) = 0;
    a1pos = 4;
    N = 1:numel(atmN);N=N';
    a = [C{a1pos}(N),C{a1pos+2}(N),C{a1pos+4}(N),C{a1pos+6}(N)];
    b = [C{a1pos+1}(N),C{a1pos+3}(N),C{a1pos+5}(N),C{a1pos+7}(N)];
    c = C{a1pos+8};
end

if nargin < 3
    charge = 0;
end

if charge == 0
    N = neutralatom(atomN);
else
    N = find(((atmN == atomN) && (charge_atoms == charge))==1);
end
a = a(N, :);
b = b(N, :);
c = c(N);
q = q(:)';
Fq = repmat(c, [1, length(q)]);
%Fq = c*ones(size(q));
q = q/(4*pi);
for j=1:numel(q)
    for i=1:4
        Fq(:,j) = Fq(:,j)+a(:,i).*exp(-b(:,i)*q(j).^2);
    end
end
Fq = Fq';