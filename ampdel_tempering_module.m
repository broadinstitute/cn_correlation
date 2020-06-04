function ampdel_tempering_module(ref_dir,perm_dir,chunk_id,iters,varargin)
%AMPDEL_TEMPERING_MODULE - executable wrapper function for ampdel_tempering()
%
%   ampdel_tempering_module(REF_DIR,PERM_DIR,CHUNK_ID,ITERS,SEED)
%
% All inputs are strings for command line usage. REF_DIR is the path to the directory 
% with input files; PERM_DIR is the path to the directory where permated output files are 
% written; CHUNK_ID is a file extension to identify the chunk; ITERS is the number of 
% iterations to do per-chunk. SEED is an optional seed to put the matlab random number 
% generator in a known state. If omitted or empty, the RNG seed is selected randomly and 
% saved in 'randseed.mat' in the OUTPUT_DIR.

% load input data from files with hardwired names in perm_dir
fprintf('loading input files\n')
load(fullfile(perm_dir,'H.mat')); % 'H' struct disruption samples
load(fullfile(perm_dir,'perm_opts.mat')); % 'perm_opts' struct for passing permutation options

% convert iteration count argument to number
if exist('iters','var') && ~isempty(iters)
    if ischar(iters) 
        iters = str2double(iters);
    end
else
    if isfield(perm_opts.Niters)
        iters = perm_opts.Niters
    else
        iters = 10;
    end
end

set_verbose_level(40); %!!! put in opts

% check for optional random number generator seed argument
if length(varargin) > 0
    % do reproducible permutations using specified seed
    seed = randseed(str2num(varargin{1}));
else
    seed = randseed();
end
fprintf('RNG seed: %d\n',seed); 
save(fullfile(perm_dir,'randseed.mat'));

% do permutations
[idx_cell,stat_finals,stats] = ampdel_tempering(H,iters,perm_opts);

% save output to files
fprintf('saving output files\n');
save(fullfile(perm_dir,['idx_cell.',chunk_id,'.mat']),'idx_cell');
save(fullfile(perm_dir,['stat_finals.',chunk_id,'.mat']),'stat_finals');
save(fullfile(perm_dir,['stats.',chunk_id,'.mat']),'stats');
fprintf('annealing_permutations module complete for chunk\n');

end % function

