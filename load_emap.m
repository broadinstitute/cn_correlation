function E = load_emap(fname,varname)
%LOAD_EMAP load event map into memory
%
%   E = load_emap(FNAME)
%
% 
    if exist('varname','var')
        e = load(fname,varname,'-mat'); %!!! does this even work? (from GISTIC load_D)
    else
        e = load(fname,'-mat');
    end

    nms = fieldnames(e);
    E = getfield(e,nms{1});
    validate_emap(E);
    
    % de-sparse the matrix
    E.dat = full(E.dat);
end
