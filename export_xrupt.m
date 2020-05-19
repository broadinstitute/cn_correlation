function export_xrupt(fname,H)
%EXPORT_XRUPT save genomic disruption structure to tab-delimited table
%
%     export_xrupt(FNAME,XRUPT)
%
% FNAME is the path to the output .txt file
% XRUPT is the genomic disruption structure input to the permutation engine
%

    % due diligence 
    assert(isfield(H,'margs'));
    assert(length(size(H.margs))) == 3);
    [Nchr,Nsample,Ntype] = size(H.margs);
    assert(isfield(H,'sdesc'));
    assert(Nsample==length(H.sdesc));
    assert(isfield(H,'chrname'));
    assert(Nchr==length(H.chrname));

    %! Ntype range?
    assert(Ntype>0);
    
    %    assert(isfield(H,'pcname'));
    %assert(isfield(H,'pcx'));
    %assert(length(H.pcx)==length(H.pcname));

    % construct sample name column
    sampcol = [H.sdesc,H.sdesc]';
    % construct alteration type column
    if Ntype == 2
        altsym = {'a';'d'};
    else
        altsym = strtrim(cellstr(num2str((1:Ntype)')))
    end
    typecol = repmat(altsym,Nsample,1);

    % construct 'pclass' permutation class column
    classlab = repmat({''},size(H.sdesc));
    for i = 1:length(H.pcx)
        classlab(H.pcx{i}) = repmat(cellstr(H.pcname{i}),size(H.pcx{i}));
    end
    classcol = [classlab,classlab]';
    % construct per-chromosome disruption data columns
    chrcols = cell(size(H.chrname));
    for c=1:Nchr
        chrlamina = squeeze(H.margs(c,:,:))';
        chrcols{c} = {H.chrname{c},chrlamina(:)};
    end
    write_filtered_tabcols(fname,[],{'sample',sampcol(:)},{'pclass',classcol(:)},...
                           {'type',typecol(:)},chrcols{:});
end
