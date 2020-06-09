function H = create_xrupt(D,margs,pcx);
%CREATE_XRUPT create genomic disruption structure from processed GISTIC inputs
%
%    XRUPT = create_xrupt(D,MARGS,PCX)
%
% D is a copy number structure
% MARGS is a 3D array (Nchromosomes x Nsamples x amp/del type) of disruption scores
% PCX defines the permutation class structure of the data
% 
% The returned XRUPT is the genomic disruption structure (permutation engine input)
%

    %! save permutation input as single package
    H = struct;
    H.sdesc = D.sdesc;
    H.margs = margs;
    H.pcx = pcx;
    H.pcname = extract_pcname(D,pcx);

%{
    H.pcname = cell(size(pcx));
    %% create class names
    % loop over permutation classes
    for i=1:length(pcx);
        % start with default "classN" name
        class_name = sprintf('class%d',i);
        % see if we can find more information in D.sis
        if isfield(D,'sis') %!!! do test outside assignment loop
            j = pcx{i}(1);
            if isfield(D.sis,'disease')
                class_name = D.sis(j).disease;
            else 
                if isfield(D.sis,'gcmtype')
                    class_name = D.sis(j).gcmtype;
                end
            end
        end
        H.pcname{i} = class_name;
    end
 %}

    % create chromosome names
    ucx = unique(D.chrn);
    H.chrname = strcat('chr',num2chromosome(ucx));
end % function

