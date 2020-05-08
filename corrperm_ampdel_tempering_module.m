function corrperm_ampdel_tempering_module(ref_dir,perm_dir,chunk_number,iters)
%ANNEALING_PERMUTATIONS_MODULE - module wrapper for ANNEALING_PERMUTATIONS
%
%   annealing_permutations_module(INPUT_DIR,OUTPUT_DIR,CHUNK_NUMBER,ITERS)
%
%   All inputs are strings for command line usage. INPUT_DIR is the path
% to the directory with input files; OUTPUT_DIR is the path to the
% directory with output files; CHUNK_NUMBER is a tring identifying the chunk; ITERS
% is the number of iterations to do per-file.

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
load([ref_dir 'margs.mat']);       % 'margs_sort' chromosome disruption values (chromosomes X samples X amp/del)
load([ref_dir 'new_samples.mat']); % 'new_samples' cell array of permutation class indices
load([ref_dir 'permute_options.mat']); % 'opts' struct for passing permutation options

set_verbose_level(40); %!!! put in opts

% do permutations
[rand_margs_cell,idx_cell,stat_finals,stats] = corrperm_ampdel_tempering(margs_sort,new_samples,iters,opts);


% save output to files
fprintf('saving output files\n')
%!save([perm_dir 'rand_margs.' chunk_number '.mat'],'rand_margs_cell')
save([perm_dir 'idx_cell.' chunk_number '.mat'],'idx_cell')
save([perm_dir 'stat_finals.' chunk_number '.mat'],'stat_finals')
save([perm_dir 'stats.' chunk_number '.mat'],'stats')
fprintf('annealing_permutations module complete\n');
