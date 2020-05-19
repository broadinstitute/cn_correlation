function corrperm_ampdel_tempering_module(ref_dir,perm_dir,chunk_id,iters)
%CORRPERM_AMPDEL_TEMPERING_MODULE - module wrapper for 'basic_tempered_annealing' algorithm
%
%   basic_tempered_annealing_module(INPUT_DIR,OUTPUT_DIR,CHUNK_ID,ITERS)
%
%   All inputs are strings for command line usage. INPUT_DIR is the path
% to the directory with input files; OUTPUT_DIR is the path to the
% directory with output files; CHUNK_ID is a string identifying the chunk; ITERS
% is the number of iterations to do per-file.

% load input data from files
fprintf('loading input files\n')
load(fullfile(ref_dir,'margs.mat'));       % 'margs_sort' chromosome disruption values (chromosomes X samples X amp/del)
load(fullfile(ref_dir,'new_samples.mat')); % 'new_samples' cell array of permutation class indices

% 'opts' struct for passing permutation options
permute_options = load(fullfile(ref_dir,'permute_options.mat')); 

% convert iteration count argument to number
if exist('iters','var') && ~isempty(iters)
    % get iteration count from 4th function argument (old way)
    if ischar(iters) 
        iters = str2double(iters);
    end
else
    % get iteration parameter from permute options (new way)
    if ~isfield(permute_options,'Niters')
        iters = permute_options.Niters;
    else
        iters = 10;
    end
end

% do permutations
[rand_margs_cell,idx_cell,stat_finals,stats] = corrperm_ampdel_tempering(margs_sort,new_samples,permute_options.Niters,permute_options.opts);

% save output to files
fprintf('saving output files\n')
%!save(fullfile(perm_dir,['rand_margs.' chunk_id '.mat'],'rand_margs_cell')
save(fullfile(perm_dir,['idx_cell.',chunk_id,'.mat']),'idx_cell');
save(fullfile(perm_dir,['stat_finals.',chunk_id,'.mat']),'stat_finals');
save(fullfile(perm_dir,['stats.',chunk_id,'.mat']),'stats');
fprintf('annealing_permutations module complete for chunk\n');
