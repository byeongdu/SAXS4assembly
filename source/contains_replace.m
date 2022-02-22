function ret = contains_replace(str1, str2)
    if exist('contains', 'file') == 2
        ret = contains(str1, str2);
    else
        if iscell(str1)
            ret = zeros(size(str1));
            for i=1:numel(ret)
                ret(i) = ~isempty(strfind(str1{i}, str2));
            end
        else
            ret = ~isempty(strfind(str1, str2));
        end
    end
