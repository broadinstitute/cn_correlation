function corrperm_ampdel_lockstep_module(ref_dir,perm_dir,cycle_number,iters)
%ANNEALING_PERMUTATIONS_MODULE - module wrapper for ANNEALING_PERMUTATIONS
%
%   annealing_permutations_module(INPUT_DIR,OUTPUT_DIR,CHUNK_NUMBER,ITERS)
%
% All inputs are strings for command line usage. INPUT_DIR is the path
% to the directory with input files; PERM_DIR is the path to the
% directory where permated output files are written; CHUNK_NUMBER is 
% an extension to identify; ITERS is the number of iterations to do per-chunk.

%randomize the column #5, the sample column

% convert iteration count argument to number
if exist('iters','var') && ~isempty(iters)
    if ischar(iters) 
        iters = str2double(iters);
    end
else
    iters = 10;
end

% load input data from files
fprintf('loading input files\n')
load(fullfile(ref_dir,'margs.mat'));       % 'margs_sort' chromosome disruption values (chromosomes X samples X amp/del)
load(fullfile(ref_dir,'new_samples.mat')); % 'new_samples' cell array of permutation class indices
load(fullfile(ref_dir,'permute_options.mat')); % 'opts' struct for passing permutation options

set_verbose_level(40); %!!! put in opts

% do permutations
[rand_margs_cell,idx_cell,stat_finals] = corrperm_ampdel_lockstep(margs_sort,new_samples,iters,opts);


% save output to files
fprintf('saving output files\n')
%!save(fullfile(perm_dir,['rand_margs.',cycle_number,'.mat']),'rand_margs_cell')
save(fullfile(perm_dir,['idx_cell.',cycle_number,'.mat']),'idx_cell')
save(fullfile(perm_dir,['stats.',chunk_number,'.mat'],'stats')
save(fullfile(perm_dir,['stat_finals.',cycle_number,'.mat']),'stat_finals')
fprintf('annealing_permutations module complete\n')


