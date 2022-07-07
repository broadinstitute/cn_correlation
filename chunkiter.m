function [table,mlt] = chunkiter(E, perm_dir, options)
    % sketch of function that iterates over file chunks
    % E is an event map structure (Estruct)
    % perm_dir path to directory of permutations
    % options are optional arguments
%
    % get sizes from observed data matrix dimensions
    [Nevents,Nsamples] = size(E.dat);
    Nchr = length(unique(E.event.chrn));

    % chromosome => event mapping
    chrns_e = cell(1,Nchr);
    for c = 1:Nchr
        chrns_e{c} = find(E.event.chrn == c);
    end

    %% process optional arguments
    if ~exist('options','var')
        options = struct;
    end
    options = impose_default_value(options,'perm_file_mask','idx_cell*');
    options = impose_default_value(options,'split_eq',false);   % set to split observed===permuted between tail and body 
    options = impose_default_value(options,'pcount',1);         % pseudo count (1 or 0) 
    options = impose_default_value(options,'analyze_lineages',false);%

    %% load disruption structure 'H' from H.mat in permutation directory
    load(fullfile(perm_dir,'H.mat'));
    assert(logical(exist('H','var')));
    assert(all(strcmp({E.sample.id},H.sdesc')));

    % initialize overall statistics
    stat = chunkstat_ExE(E.dat);

    % lineage analysis based on permutation lineages
    Nlineages = 0;
    if options.analyze_lineages
        lindices = H.pcx;
        linnames = H.pcname;
        Nlineages = length(H.pcx);
    end

    %% initialize statistics

    % initiate per-lineage statistics
    if Nlineages > 0
        lstats = cell(1,Nlineages);
        for l = 1:Nlineages
            lstats{l} = chunkstat_ExE(E.dat(:,lindices{l}));
        end
    end

    
    %% process permutation directory
    files = dir(fullfile(perm_dir,options.perm_file_mask));
    
    Nfiles = length(files);
    verbose('Reading %d permutation chunks from ''%s''',10,Nfiles,perm_dir);
    if ~Nfiles
        error('No permutation data to process');
    end

    %% loop over files
    Nperms = 0;
    for k = 1:Nfiles
        verbose('Processing ''%s''',20,files(k).name);
        tic
        load(fullfile(perm_dir,files(k).name));  % 'idx_cell' cell array, each element Nchr x Nsamples
        if ~exist('NExpectedPerms','var')
            npf = length(idx_cell);
            NExpectedPerms = Nfiles * npf;
            hash32 = zeros(NExpectedPerms,1);
        end

        % loop over permutations in the chunk
        for i = 1:npf
            % count permutation
            Nperms = Nperms + 1;        

            % use permutation indices to rearrange observed events
            idx_mat = idx_cell{i};
            % check for old 3D index matrix and convert if necessary
            if length(size(idx_mat))==3
                idx_mat = idx_mat(:,:,1);
            end

            % calculate CRC as test for identical permutations
            hash32(Nperms) = crc32(uint32(idx_mat));

            % map chromosomes to events one chromosome at a time
            emap = zeros(Nevents,Nsamples); % allocate storage for expanded event map
            for c = 1:Nchr
                emap(chrns_e{c},:) = double(E.dat(chrns_e{c},idx_mat(c,:)));
            end
            
            % update statistics with current permutation
            stat = perm(stat,emap);
            for l = 1:Nlineages
                lstats{l} = perm(lstats{l},emap(:,lindices{l}));
            end
            
        end % loop over permutations in chunk
    end % loop over chunk files
    
    % calculate number of pairs on different chromosomes
    Npairs = 0;
    for i = 1:Nevents
        Npairs = Npairs + sum(E.event.chrn(i) < E.event.chrn);
    end
    % create indexed pair of indices
    pairs_idx = zeros(Npairs,2);
    s = 1;
    for i = 1:Nevents-1
        for j = i+1:Nevents
            if E.event.chrn(i) ~= E.event.chrn(j)
                pairs_idx(s,:) = [i,j];
                s = s + 1;
            end
        end
    end

    if length(unique(hash32)) ~= Nperms
        warning('some permutations might be identical');
        verbose('%d of %d hashes unique',10,length(unique(hash32)),Nperms);
    end
    
    %% results
    
    % pair name columns
    e1 = E.event.name(pairs_idx(:,1));
    e2 = E.event.name(pairs_idx(:,2));
    %! paircols = struct('event1',e1,'event2',e2);

    % overall
    table = results(stat,pairs_idx,options);
    [table.event1] = deal(e1{:});
    [table.event2] = deal(e2{:});
    table = orderfields(table,{'event1','event2','p_corr','p_anti'})

    % melted lineage table
    ltabs = cell(Nlineages,1);
    for l = 1:Nlineages
        tab = results(lstats{l},pairs_idx,options);
        [tab.lineage] = deal(linnames{l});
        [tab.event1] = deal(e1{:});
        [tab.event2] = deal(e2{:});
        ltabs{l} = tab;
    end
    mlt = vertcat(ltabs{:});
    mlt = orderfields(mlt,{'lineage','event1','event2','p_corr','p_anti'})
        
end % function
