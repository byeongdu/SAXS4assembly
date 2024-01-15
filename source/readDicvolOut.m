function cellp = readDicvolOut(filename)
% this function parses the output of DICVOL software.
alllines = readlines(filename);
ltype = 0;
A = NaN;
B = NaN;
C = NaN;
ALPHA = 90;
BETA = 90;
GAMMA = 90;
structures = {};
strc = [];
str_index = 0;
factor = 1;
toskip = false;
t_format = false;
for i=1:numel(alllines)
    l = alllines{i};
    if contains(l, 'with factor')
        dt = split(l, "with factor");
        dt = split(strtrim(dt{2}), ' ');
        factor = str2double(dt{1});
    end
    if contains(l, 'C U B I C    S Y S T E M')
        ltype = 'cubic';
    end
    if contains(l, 'T E T R A G O N A L   S Y S T E M')
        ltype = 'tetragonal';
    end
    if contains(l, 'H E X A G O N A L    S Y S T E M')
        ltype = 'hexagonal';
    end
    if contains(l, 'O R T H O R H O M B I C    S Y S T E M')
        ltype = 'orthorhombic';
    end
    if contains(l, 'M O N O C L I N I C    S Y S T E M')
        ltype = 'monoclinic';
    end
    if contains(l, 'T R I C L I N I C    S Y S T E M')
        ltype = 'triclinic';
    end
    if contains(l, 'THE SOLUTION IS NOW USED TO TRY TO INDEX ALL INPUT')
        toskip = true;
    end
    if t_format
        dt = split(l, "VOL=");
        vol = str2double(dt{2});
        p_str = strtrim(dt{1});
        strv = replace(p_str, ' ', ';');
        strv = replace(strv, ';;', ';');
        strv = replace(strv, '=;', '=');
        eval([strv, ';'])
        ALPHA = ALP;
        BETA = BET;
        GAMMA = GAM;
        t_format = false;
        strc.latticetype = ltype;
        strc.cellp = [A/factor, B/factor, C/factor, ALPHA, BETA, GAMMA];
        if str_index == 0
            struct0 = strc;
        end
        str_index = str_index + 1;
        structures{str_index} = strc;
    end
    if contains(l, 'DIRECT PARAMETERS :') 
        if toskip
            toskip = false;
            continue;
        end
        dt = split(l, "DIRECT PARAMETERS :");
        dt = split(dt{2}, "VOLUME=");
        vol = str2double(dt{2});
        p_str = strtrim(dt{1});
        strv = split(p_str, '  ');
        for k=1:2:numel(strv)
            eval([strv{k}, strv{k+1}, ';']);
        end
        % if contains(p_str, "A=")
        %     strv = split(p_str, "A=");
        %     cell_index = 1;
        % end
        % 
        switch ltype
            case 'cubic'
                B = A;
                C = A;
                ALPHA = 90;
                BETA = 90;
                GAMMA = 90;
            case 'hexagonal'
                B = A;
                ALPHA = 90;
                BETA = 90;
                GAMMA = 120;
            case 'tetragonal'
                B = A;
                ALPHA = 90;
                BETA = 90;
                GAMMA = 90;                
            case 'orthorhombic'
                ALPHA = 90;
                BETA = 90;
                GAMMA = 90;
            case 'monoclinic'
                ALPHA = 90;
                GAMMA = 90;
        end
    
        strc.latticetype = ltype;
        strc.cellp = [A/factor, B/factor, C/factor, ALPHA, BETA, GAMMA];
        if str_index == 0
            struct0 = strc;
        end
        str_index = str_index + 1;
        structures{str_index} = strc;
    end
    if contains(l, 'DIRECT PARAMETERS AND THEIR STANDARD DEVIATIONS :')
        t_format = true;
    end
end
if str_index == 1
    cellp = struct0;
else
    cellp = structures;
end