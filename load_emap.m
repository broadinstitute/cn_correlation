function E = load_emap(fname,varname)
%LOAD_EMAP load event map into memory
%
%   E = load_emap(FNAME)
%
% 
    if exist('varname','var')
        e = load(fname,varname); %!!! does this work? (from GISTIC load_D)
    else
        e = load(fname);
    end

    nms = fieldnames(e);
    E = getfield(e,nms{1});
    % validate
    assert(isfield(E,'sdesc'));
    assert(isfield(E,'dat'));
    assert(isfield(E,'event'));
    assert(isfield(E.event,'name'));
    assert(isfield(E.event,'type'));
    assert(isfield(E.event,'chrn'));
    %@!!! validate sizes

    % de-sparse the matrix
    E.dat = full(E.dat);
end
