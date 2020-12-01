function E = validate_emap(E)
    assert(isfield(E,'sample'));
    assert(isfield(E.sample,'id'));
    assert(isfield(E,'dat'));
    assert(isfield(E,'event'));
    assert(isfield(E.event,'name'));
    assert(isfield(E.event,'type'));
    assert(isfield(E.event,'chrn'));
    %@!!! validate sizes
end
