function save_emap(E,fname)
%SAVE_EMAP save event map in matlab native format
%
%    save_emap(E,FNAME)
%

    % validate
    validate_emap(E);

    % compress sparse enough call data
    if sum(sum(E.dat)) / (size(E.dat,1) * size(E.dat,2)) < 0.06
        E.dat = sparse(E.dat);
    end
    save(fname,'E');
end % function
