function I = FileIter2(perm_dir)
%% create file iterator by reading perm_opts.mat from permutation directory
%
    
% closure data
    chunk = [];  % current chunk
    iter = [];   % position in current chunk
    idx_cell = []; % chunk data (cell array of permutation indices)
    perm_path_template = []; % sprintf string for chunk file path
    chunkx = []; % index to chunks backed by actual files
    NactualChunks = []; % number of chunks backed by actual files
    Npermutations = []; % total number of permutations
    
    permcount = []; % interator index
    
    I = struct;
    I.options = [];

    x = load(fullfile(perm_dir,'perm_opts.mat')); % load perm_opts
    I.options = x.perm_opts;
    % create a map of permutations that actually exist
    done_map = false(1,I.options.Nchunks);
    perm_path_template = fullfile(perm_dir,I.options.output_subdir,I.options.chunk_template);
    % find chunks whose outputs have begun being written
    for i = 1:I.options.Nchunks
        done_map(i) = logical(exist(sprintf(perm_path_template,i),'file'));
    end
    chunkx = find(done_map);
    NactualChunks = length(chunkx);
    Npermutations = NactualChunks*I.options.Niters;
    %!!! done_map and chunkx are not currently used (existence test at chunk load instead)
    
    reset();

    function n = count()
        n = Npermutations;
    end
        
    %% restart iterations at beginning
    % no status refresh from file system
    function reset()
        chunk = 0;
        iter = I.options.Niters + 1; %!!! hacky
        idx_cell = [];
        permcount = 1;
    end

    %% pull next permutation
    function perm_idx = next()
    %
    % load chunk if needed
        if (iter > I.options.Niters) && (chunk < I.options.Nchunks)
            load_next_chunk();
        end
        if iter <= I.options.Niters
            perm_idx = idx_cell{iter};
            iter = iter + 1;
            permcount = permcount + 1;
        else
            perm_idx = [];
        end
    end

    %% are there any more permutations?
    function tf = any_more()
        tf = permcount <= Npermutations;
%!        tf = (chunk < I.options.Nchunks) || (iter <= I.options.Niters);
    end
    
    %!    end % methods

    %! methods(Access = 'private')
    %% load a chunk of permutations
    function load_next_chunk()
        chunk = chunk + 1;
        file = sprintf(perm_path_template,chunk);
        while ~exist(file,'file') && chunk <= I.options.Nchunks
            verbose('skipping missing chunk %d',10)
            chunk = chunk + 1
        end
        if chunk <= I.options.Nchunks
            verbose('loading chunk %d',20,chunk);
            inbox = load(file);
        end
        idx_cell = inbox.idx_cell;
        iter = 1;
    end

    % return function handles
    I.count = @count;
    I.reset = @reset;
    I.next = @next;
%!  I.load_chunk = @load_chunk;
    I.any_more = @any_more;
end

