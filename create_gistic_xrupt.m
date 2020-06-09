function H = create_gistic_xrupt(D,margs,pcx);
%CREATE_GISTIC_XRUPT create genomic disruption structure from processed GISTIC inputs
%
%    XRUPT = create_gistic_xrupt(D,MARGS,PCX)
%
% inputs are outputs of corrperm_prep:
%   D is a copy number structure
%   MARGS is a 3D array (Nchromosomes x Nsamples x amp/del type) of disruption scores
%   PCX defines the permutation class structure of the data
% 
% The returned XRUPT is the genomic disruption structure (permutation engine input)
%
    %! save permutation input as single package
    H = struct;
    H.sdesc = D.sdesc;
    H.margs = margs;
    H.pcx = pcx;
    H.pcname = extract_pcname(D,pcx);

    % create chromosome names
    ucx = unique(D.chrn);
    H.chrname = strcat('chr',num2chromosome(ucx));
end % function

