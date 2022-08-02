function x = perm_status(perm_dir)
% prototype permutation status function    
    load(fullfile(perm_dir,'perm_opts.mat')); % load perm_opts

    Nchunks = perm_opts.Nchunks;
    perm_output_dir = fullfile(perm_dir,perm_opts.output_subdir);

    % initialize chunk maps
    done_map = false(1,Nchunks);
    start_map = false(1,Nchunks);
    err_map = false(1,Nchunks);
    iter_count = zeros(1,Nchunks);
    
    for i = 1:Nchunks
        % find chunks whose outputs have begun being written
        if exist(fullfile(perm_output_dir,sprintf(perm_opts.chunk_template,i)),'file')
            done_map(i) = true;
        end
        % find chunks that have started running enough to generate output (may not work on all MPEs)
        file_info = dir(fullfile(perm_output_dir,sprintf(perm_opts.stdout_template,i)));
        outfile = fullfile(perm_output_dir,file_info.name);
        if ~isempty(file_info) && file_info.bytes > 0
            start_map(i) = true;
            [err,out] = unix(sprintf('cat %s | grep -E "iteration [0-9]+:" | wc -l',outfile));
            if err==0
                iter_count(i) = str2num(strtrim(out));
            end
        end
        % find chunks with errors
        errfile = dir(fullfile(perm_output_dir,sprintf(perm_opts.stderr_template,i)));
        if ~isempty(errfile) && errfile.bytes > 0
            err_map(i) = true;
        end
    end
    % test output
    x = struct;
    x.done_map = done_map;
    x.err_map = err_map;
    x.start_map = start_map;
    x.iter_count = iter_count;
    
end
