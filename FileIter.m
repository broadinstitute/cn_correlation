% implementation of a permutation iterator for chunk files
classdef FileIter
    properties
        % set at construction
        perm_path_template = '';
        Nchunks = 0;
        Niters = 1;
        
        % dynamic
        chunk = 0; % chunk index
        iter = 0; % permutation within chunk index
        idx_cell = [];
    end

    methods
        %% constructor called before chunk processing
        % I = FileIter(perm_dir,file_opts))
        % emap is an Nevents x Nsamples logical array of observed events in each sample
        function I = FileIter(perm_dir)
            load(fullfile(perm_dir,'perm_opts.mat')); % load perm_opts
            file_opts = perm_opts;
            I.Nchunks = file_opts.Nchunks;
            I.Niters = file_opts.Niters;
            I.perm_path_template = fullfile(perm_dir,file_opts.chunk_template);
            I = reset(I);
        end
        
        %% restart iterations at beginning
        % no status refresh from file system
        function I = reset(I)
            I.chunk = 0;
            I.iter = I.Niters + 1; %!!! hacky
            I.idx_cell = [];
        end

        %% pull next permutation
        function perm_idx = next(I)
        %
            % load chunk if needed
            if I.iter > I.Niters && I.chunk < I.Nchunks
                load_chunk(I,I.chunk+1)
            end
            if I.iter <= I.Niters
                perm_idx = I.idx_cell{I.iter};
                I.iter = I.iter + 1;
            else
                perm_idx = [];
            end
        end

        %% are there any more permutations?
        function tf = any_more(I)
            tf = I.iter <= I.Niters && I.chunk <= I.Nchunks;
        end
        
    end % methods

    methods(Access = 'private')
        %% load a chunk into iterator
        function load_chunk(I,k)
            inbox = load(sprintf(I.perm_path_template,k));
            I.idx_cell = inbox.idx_cell;
            I.chunk = k;
            I.iter = 1;
        end        
    end
end
