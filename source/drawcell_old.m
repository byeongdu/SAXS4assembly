function [App, Ap, mat, phandle] = drawcell(varargin)

%function [App, Ap, mat] = drawcell(sginfo, atom, A, B, C, AL, BE, GA, numofcell, nodisplay)
% function [App, Ap, mat] = drawcell(sginfo, atom, A, B, C, AL, BE, GA,
% numofcell, nodisplay)
% function [App, Ap, mat] = drawcell(sginfo, atom, cellinfo,numofcell, nodisplay)
% function [App, Ap, mat] = drawcell(sginfo, atom, cellinfo,numofcell)
% function [App, Ap, mat] = drawcell(sginfo, atom, cellinfo)
%
% if a particle is composed of atoms, use drawcell_atm.m
% atom = [1 x y z; 2 x y z; 2 x y z; 3 x y z;...]
% example :
%   an = drawcell(sginfo, [0.3, 0, 0, 0;0.2, 0.3333, 0.6667, 0;0.2, 0.5, 0, 0.5],
%   1, 1, 0.8029, 90, 90, 120);
%
% Format for atom or SgAtoms
% [radius, positionx, positiony, positionz, SOF, TfactorB, Lshell,
% rho_core, rho_shell, rho_solv];
% 
latticevectors = [];
roundres = 1E15;
phandle = [];
numofcell = 1;
nodisplay = 0;
sginfo = varargin{1};
atom = varargin{2};
if numel(varargin) == 10
    sginfo = varargin{1};
    atom = varargin{2};
    A = varargin{3};
    B = varargin{4};
    C = varargin{5};
    AL = varargin{6};
    BE = varargin{7};
    GA = varargin{8};
    numofcell = varargin{9};
    nodisplay = varargin{10};
end    
if numel(varargin) == 9
    sginfo = varargin{1};
    atom = varargin{2};
    A = varargin{3};
    B = varargin{4};
    C = varargin{5};
    AL = varargin{6};
    BE = varargin{7};
    GA = varargin{8};
    numofcell = varargin{9};
end    
if numel(varargin) == 8
    sginfo = varargin{1};
    atom = varargin{2};
    A = varargin{3};
    B = varargin{4};
    C = varargin{5};
    AL = varargin{6};
    BE = varargin{7};
    GA = varargin{8};
end    
if numel(varargin) < 6
    sginfo = varargin{1};
    atom = varargin{2};
    if numel(varargin) > 2
        cellinfo = varargin{3};
        A = cellinfo.A;
        B = cellinfo.B;
        C = cellinfo.C;
        AL = cellinfo.alpha;
        BE = cellinfo.beta;
        GA = cellinfo.gamma;
        latticevectors = cellinfo.latticevectors;
        if numel(varargin)>3
            numofcell = varargin{4};
        end
        if numel(varargin)>4
            nodisplay = varargin{5};
        end
    end
end    
if ~isstruct(sginfo)
    if ischar(sginfo)
        sginfo = sgroup(sginfo);
    end
end
if ~isfield(sginfo, 'nsymop_prim')
    % this new spacegroup has all symoperation including lattice centering
    % and inversion. Let's use this. For now, will convert old one to new..
    sginfo = sg_gen(sginfo.Number);
end

NosymMtx = sginfo.NoSymMatrices;

if iscell(atom)
    NoAtom = numel(atom);
    pt = atom;
else
    list_atoms = atom(:,1);

    [LAtom, ord] = unique(atom(:,1));
    LAtom = flipud(LAtom);
    ord = flipud(ord);

    cln = numel(atom(1,:));
    %Numinput = numel(atom(:,1));
    k = 0;
    for i=1:numel(LAtom)
        k = k+1;
        inx = find(list_atoms == LAtom(i));

        TfactorB = 0;SOF = [];Lshell = [];rho_shell = [];rho_solv = [];rho = 1;
        pt{k}.radius = atom(inx(1),1);
        pt{k}.position = atom(inx(1),2:4);
        pt{k}.index = inx(1);
        if cln > 4
            SOF = atom(pt{k}.index, 5);
        end
        if cln > 5
            TfactorB = atom(pt{k}.index, 6);
        end
        if cln > 6
            Lshell = atom(pt{k}.index, 7);
            rho_core = atom(pt{k}.index, 8);
            rho_shell = atom(pt{k}.index, 9);
            rho_solv = atom(pt{k}.index, 10);
            rho = [rho_core, rho_shell, rho_solv];
        end
        pt{k}.TfactorB = TfactorB;
        pt{k}.SOF = SOF;
        pt{k}.Lshell = Lshell;
        pt{k}.rho = rho;

        if numel(inx) > 1
            for m=2:numel(inx)
                nw = atom(inx(m), 5:end);
                isitdifferent = 0;
                for mm=1:m-1
                    tsub = find(atom(inx(mm), 5:end)-nw ~= 0, 1);
                    if ~isempty(tsub)
                        isitdifferent = 1;
                    end
                end
                if isitdifferent == 1
                    k = k+1;
                    TfactorB=0;SOF = [];Lshell = [];rho_shell = [];rho_solv = [];rho = 1;
                    pt{k}.radius = atom(inx(m), 1);
                    pt{k}.position = atom(inx(m), 2:4);
                    pt{k}.index = inx(m);
                    if cln > 4
                        SOF = atom(pt{k}.index, 5);
                    end
                    if cln > 5
                        TfactorB = atom(pt{k}.index, 6);
                    end
                    if cln > 6
                        Lshell = atom(pt{k}.index, 7);
                        rho_core = atom(pt{k}.index, 8);
                        rho_shell = atom(pt{k}.index, 9);
                        rho_solv = atom(pt{k}.index, 10);
                        rho = [rho_core, rho_shell, rho_solv];
                    end
                    pt{k}.TfactorB = TfactorB;
                    pt{k}.SOF = SOF;
                    pt{k}.Lshell = Lshell;
                    pt{k}.rho = rho;
                else
                    pt{k}.position = [pt{k}.position;atom(inx(m), 2:4)];
                end
            end
        end
    end
    NoAtom = k;
end

% cell coordinate calculation.

for k=1:NoAtom
    Ap{k} = pt{k};
    Ap{k}.position = [];
    pos = pt{k}.position;
    num_pos = numel(pos)/3;
    for i=1:sginfo.
        R = sginfo.symop(i).rot;
        mat = Hmat*R*Hmat2;
        % mat : matrix for vertices in Cartesian coordinate.
        T = sginfo.symop(i).trn; T = T(:);
        % 1. Transform vertices..
        vertn = R*vert + T;
    
    for pa = 1:num_pos
        tt = [];bt = [];
        pos_curr = pos(pa, :);
        for i=1:NosymMtx
            tt(i,:) = pos_curr*sginfo.SymMatrices(:,1:3,i)' + sginfo.SymMatrices(:,4,i)'/sginfo.Tvect_scale;
            if sginfo.Flag_1bar == 1
                bt(i,:) = -pos_curr*sginfo.SymMatrices(:,1:3,i)' + sginfo.SymMatrices(:,4,i)'/sginfo.Tvect_scale;
            end
        end
        tt = [tt;bt];
        %tt = reducehkl(round(tt*roundres))/roundres;
        [~, ind] = reducehkl(tt);
        tt = tt(ind, :);
        t = find(tt < 0); tt(t) = tt(t) + 1;
        t = find(tt >= 1); tt(t) = tt(t) - 1;
        t = find((tt(:,1) <= -1) | (tt(:,2) <= -1) | (tt(:,3) <= -1));
        tt(t,:) = [];
        t = find((tt(:,1) >= 1) | (tt(:,2) >= 1) | (tt(:,3) >= 1));
        tt(t,:) = [];
        % lattice centering 
        % P, F, I, B
        atmm = tt;
        Natm = numel(atmm)/3;
        if size(sginfo.SymMatrices, 3)<24
            NoLatticeCenteringVector = 1; % R centering;
        else
            NoLatticeCenteringVector = sginfo.NoLatticeCenteringVector;
        end
        for lc = 1:NoLatticeCenteringVector
            for i=1:Natm;
                tt(i + (lc-1)*Natm,:) = atmm(i,:) + sginfo.LatticeCenteringVector(:,lc)';
            end;
        end
        %tt=reducehkl(round(tt*roundres))/roundres;
        [~, ind] = reducehkl(tt);
        tt = tt(ind, :);
        
        t = find(tt < 0); tt(t) = tt(t) + 1;
        t = find(tt >= 1); tt(t) = tt(t) - 1;
        t = find((tt(:,1) <= -1) | (tt(:,2) <= -1) | (tt(:,3) <= -1));
        tt(t,:) = [];
        t = find((tt(:,1) >= 1) | (tt(:,2) >= 1) | (tt(:,3) >= 1));
        tt(t,:) = [];
        Ap{k}.position=[Ap{k}.position;tt];
    end
%    Ap{k}.position = reducehkl(round(Ap{k}.position*roundres))/roundres;
    [~, ind] = reducehkl(Ap{k}.position);
    Ap{k}.position = Ap{k}.position(ind, :);
end

if nargin == 2;
    App = Ap;
    return;
end

% Translation
for k=1:NoAtom
    atm = Ap{k}.position;
    nc = 1;
    %atm = [atm; ...
    %    [atm(:,1)+nc, atm(:,2), atm(:,3)];...
    %    [atm(:,1), atm(:,2)+nc, atm(:,3)];...
    %    [atm(:,1), atm(:,2), atm(:,3)+nc];...
    %    [atm(:,1)+nc, atm(:,2)+nc, atm(:,3)];...
    %    [atm(:,1)+nc, atm(:,2), atm(:,3)+nc];...
    %    [atm(:,1), atm(:,2)+nc, atm(:,3)+nc];...
    %    [atm(:,1)+nc, atm(:,2)+nc, atm(:,3)+nc];...
    %    ];
    natmall = [];
    if numel(numofcell) == 1
        if A > 0
            numofcellx = numofcell;
        else
            numofcellx = 0;
        end
        if B > 0
            numofcelly = numofcell;
        else
            numofcelly = 0;
        end
        if C > 0
            numofcellz = numofcell;
        else
            numofcellz = 0;
        end
    else
        if A > 0
            numofcellx = numofcell(1);
        else
            numofcellx = 0;
        end
        if B > 0
            numofcelly = numofcell(2);
        else
            numofcelly = 0;
        end
        if C > 0
            numofcellz = numofcell(3);
        else
            numofcellz = 0;
        end
    end

    for nc1 = -fix(numofcellx/2):ceil(numofcellx/2)
        for nc2 = -fix(numofcelly/2):ceil(numofcelly/2)
            for nc3 = -fix(numofcellz/2):ceil(numofcellz/2)
                natm = [atm(:,1)+nc1, atm(:,2)+nc2, atm(:,3)+nc3];
                natmall = [natmall;natm];
            end
        end
    end
    atm = natmall;
    t = find((atm(:,1) <-fix(numofcellx/2)) | (atm(:,2) <-fix(numofcelly/2)) | (atm(:,3) <-fix(numofcellz/2)));
    atm(t,:) = [];
    t = find((atm(:,1) > (ceil(numofcellx/2))) | (atm(:,2) > (ceil(numofcelly/2))) | (atm(:,3) > (ceil(numofcellz/2))));
    atm(t,:) = [];
    %atm = reducehkl(round(atm);
    %atm = reducehkl(round(natmall);
    Ap{k}.position = atm;
end
rd = pi/180;
% laboratory realspace coordinate..
sginfo.LatticeSystem = strtrim(sginfo.LatticeSystem);
if isempty(latticevectors)
    switch sginfo.LatticeSystem
        case 'trigonal'
            ax = [A, 0, 0];
            bx = [A*cos(pi/180*120), A*sin(pi/180*120), 0];
            cx = [0, 0, C];        
        case 'hexagonal'
            ax = [A, 0, 0];
            bx = [A*cos(pi/180*120), A*sin(pi/180*120), 0];
            cx = [0, 0, C];        
        case 'cubic'
            ax = [A, 0, 0];
            bx = [0, A, 0];
            cx = [0, 0, A];
        case 'tetragonal'
            ax = [A, 0, 0];
            bx = [0, B, 0];
            cx = [0, 0, C];
        case 'orthorhombic'
            ax = [A, 0, 0];
            bx = [0, B, 0];
            cx = [0, 0, C];
        otherwise
            ax = [A, 0, 0];
            bx = [B*cos(rd*GA), B*sin(rd*GA), 0];
            tmp = (cos(rd*AL)-cos(rd*GA)*cos(rd*BE))/sin(rd*GA);
            tmp2 = sqrt(sin(rd*BE)^2*sin(rd*GA)^2-(cos(rd*AL)-cos(rd*GA)*cos(rd*BE)));
            cx = [C*cos(rd*BE), C*tmp, C*tmp2];
    end
else
    ax = latticevectors(1, :);
    bx = latticevectors(2, :);
    cx = latticevectors(3, :);
end
mat = [ax',bx',cx'];

for k=1:NoAtom
    %App{k}.position = reducehkl(round(round(Ap{k}.position*mat'*1000)/1000);
    App{k}.position = Ap{k}.position*mat';
end
%figure;
if nodisplay == 1
    return
end
[x, y, z] = sphere(20);
for k=1:NoAtom
    if isfield(pt{k}, 'color')
        cl = pt{k}.color;
    else
        switch k
            case 1
                cl = 'g';
            case 2
                cl = 'b';
            case 3
                cl = 'r';
            case 4
                cl = 'k';
            case 5
                cl = 'm';
            otherwise
                cl = 'y';
        end
    end
    if isfield(pt{k}, 'shape') 
        switch pt{k}.shape
            case 'sphere'
                %pt{k}.radius = pt{k}.edgelength;
            otherwise
                pt{k}.radius = pt{k}.edgelength/2;
        end
    end
           
            
    %visualsphere(LAtom(k),App{k}.position, cl, 20);hold on; % draw core.
    ph = visualsphere(pt{k}.radius, App{k}.position, cl, 20);hold on; % draw core.
    pshandle = [];
    App{k}.position = reducehkl(App{k}.position);
    if ~isempty(pt{k}.Lshell)
        Lsh = pt{k}.Lshell+pt{k}.radius;
        pshandle = visualsphere(Lsh,App{k}.position, cl, 20, 1); % draw shells if they are available
        %visualsphere(Lshell(k),App{k}.position, cl, 20, 1); % draw shells if they are available
    end
    ph = [ph, pshandle];
    
    % set userdata.
    for t = 1:numel(App{k}.position)/3
        pinfo.position = App{k}.position(t, :);
        pinfo.positionincell = Ap{k}.position(t, :);
        pinfo.particle = pt{k};
        pinfo.wk = pt{k}.position;
        set(ph(t), 'userdata', pinfo)
    end
    phandle = [phandle, ph];
end
% lighting phong; camlight left;
% draw unit cell
% cellc = [0,0,0;1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;1,1,1];
% ncellc = [cellc(:,1)*A, cellc(:,2)*B, cellc(:,3)*C];
%drawunitcell(cellinfo)
% ncellc = cellc*mat';
% for k=1:8
%     for j=k+1:8
%         dn = norm(cellc(k,:)-cellc(j,:));
%         if dn == 1
%             line([ncellc(k,1), ncellc(j,1)], [ncellc(k,2), ncellc(j,2)],[ncellc(k,3), ncellc(j,3)])
%         end
%     end
% end
% xlabel('a')
% ylabel('b')
% zlabel('c')

%clf;visualsphere(0.1, t1, 'r',20)