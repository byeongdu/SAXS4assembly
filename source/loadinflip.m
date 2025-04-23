function [cellinfo, sg, fn, celldm, data, fit_info] = loadinflip(varargin)
    path_ = fileparts(varargin{1});
    ff = fopen(varargin{1});
    shfactor = 1;
    axv = [];
    bxv = [];
    cxv = [];
    data = [];
    fit_info = [];
    rd = 0;
    while (1)
        fL = fgetl(ff);
        if ~ischar(fL), break, end
        [aa, bb] = strtok(fL, ' ');
        switch aa
            case 'outputfile'
                if isempty(path_)
                    fn = strtrim(bb);
                else
                    fn = [path_, filesep, strtrim(bb)];
                end
            case 'cell'
                cellp = eval(strcat('[',strtrim(bb),']'));
                if numel(cellp)>3 % for 3D
                    cellp(1:3) = cellp(1:3)*shfactor; %convert A into nm.
                    celldm = 3;
                elseif numel(cellp)==3 % for 2D
                    %cell(1:2) = cell(1:2)*shfactor/10; %convert A into nm.
                    cellp = [cellp(1:2)*shfactor, 0, 90, 90, cellp(3)];
                    celldm = 2;
                end
            case '#latticevector(a)'
                axv = eval(strcat('[',strtrim(bb),']'));
            case '#latticevector(b)'
                bxv = eval(strcat('[',strtrim(bb),']'));
            case '#latticevector(c)'
                cxv = eval(strcat('[',strtrim(bb),']'));
            case '#cellshrinkingfactor'
                shfactor = str2double(strtrim(bb));
            case '#Gaussian_Width'
                fit_info = setfield(fit_info, 'Gaussian_Width', str2double(bb));
            case '#Domain_Size'
                fit_info = setfield(fit_info, 'Domain_Size', str2double(bb));
            case '#Microstrain'
                fit_info = setfield(fit_info, 'Microstrain', str2double(bb));
            case '#SF_userBG'
                fit_info = setfield(fit_info, 'SF_userBG', str2double(bb));                
            case '#spacegroup'
                sgstr = strtrim(bb);
                if ~isempty(sgstr)
                    %sg = sgroup(sgstr);
                    sg = sg_gen(sgstr);
                    wy = wyckoff;
                    if isfield(sg, 'Number')
                        Nsg = sg.Number;
                    end
                    if isfield(sg, 'number')
                        Nsg = sg.number;
                    end
                    if Nsg > 0
                        sg.wyckoff = wy{Nsg};
                    end
                    assignin('base', 'ccp4sg', sg)
                end
            case 'fbegin'
                rd = 1;
            case 'endf'
                rd = 0;
        end

        if (rd == 1) & ~strcmp(aa,'fbegin') & ~strcmp(aa,'endf')
            if strcmp(fL(1), '#')
                continue
            end
            data = [data; eval(strcat('[',fL, ']'))];
        end                
    end
    if ~isempty(axv)
        cellinfo = celcon4latvec(axv, bxv, cxv);
    else
        cellinfo = celcon(cellp);
    end
    assignin('base', 'cellinfo', cellinfo)

    fclose(ff);
end
