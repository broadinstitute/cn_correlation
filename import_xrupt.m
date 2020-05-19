function H = import_xrupt(fname)
%IMPORT_XRUPT read genomic disruption structure from text file
%
%    XRUPT = import_xrupt(FNAME)
%
% FNAME is the path to the output .txt file
% XRUPT is the genomic disruption structure input to the permutation engine
%

    % read in columns using R table utility
    dtab = read_R_table(fname);
    assert(isfield(dtab,'sample'));
    assert(isfield(dtab,'type'));
    assert(isfield(dtab,'pclass'));
    % read in samples (duplicates because each type has own line)
    [sdesc,id,ix] = unique({dtab.sample},'stable');
    Nsample = length(sdesc);
    % read in types
    [upc,~,i2] = unique({dtab.type},'stable');
    Ntype = length(upc);

    % make sure we have every (sample,type) element exactly once
    aa = accumarray([ix,i2],1);
    assert(all(all(aa==1))); %@!!! missing data error

    % create sample extraction index for pulling Nsample x Ntype lamina out of each column
    xtrax = accumarray([ix,i2],(1:(Nsample*Ntype)'));

    % get chromosome names from field names using field order
    chrname = setdiff(fieldnames(dtab),{'sample','type','pclass'},'stable');
    Nchr = length(chrname);
    % build 3D disruption matrix 
    margs = zeros([Nchr,Nsample,Ntype]);
    % loop over chromosomes
    for c=1:Nchr
        % get vector of chromosome data
        chrdat = [dtab.(chrname{c})];
        % sort it into Nsample x Ntype lamina disruption matrix
        margs(c,:,:) = chrdat(xtrax);
    end

    % permutation class processing
    [pcname,~,ip] = unique({dtab(id).pclass});
    Nclass = length(pcname);
    %! NOTE: permutation class comes from first type of each sample
    %! could enforce consistency across types
    pcx = cell(1,Nclass);
    for clx = 1:Nclass
        pcx{clx} = find(clx==ip)';
    end

    % build output structure
    H = struct;
    H.sdesc = sdesc';
    H.margs = margs;
    H.pcx = pcx;
    H.pcname = pcname;
    H.chrname = chrname';
end

