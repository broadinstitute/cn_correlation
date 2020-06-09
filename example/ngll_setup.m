% script to define global variables for "low level" analysis of 2013 Nature Genetics pan-cancer study

clear

ext = 'ngll'; % dinstiguishing file extension
ref_dir = fullfile(pwd,'ngll_work/');     % where permutation inputs are stored
perm_dir = fullfile(ref_dir,'ll_permout'); % where permutation output chunk files are stored
results_dir = fullfile(ref_dir,'/ll_results/'); % where significance results are stored

%% load input data for permutations and analysis
% load disruption profile for permutations
H = import_xrupt('ngll.xrupt.txt');
% load event map for analysis
%E = import_emap('ngll.emap.txt'); %!!! fix emap file format
E = load_emap('ngll.emap.mat');

%% set parameters for low-level tempered annealing
perm_opts = struct;
% misc
perm_opts.algorithm = 'basic_tempered_annealing'; % (for validation, does not set algorithm)
perm_opts.stat_interval = 1000;
perm_opts.class_subiters = 100;
perm_opts.ramp_steps = 10;
perm_opts.minbinmin = 0;
% initial cycle
perm_opts.initial.iters = 3e4;
perm_opts.initial.step = 10;
perm_opts.initial.temp = [1e5,1e4];
% single tempering cycle
perm_opts.cycle(1).iters = 5000;
perm_opts.cycle(1).step = 10;
perm_opts.cycle(1).heatup = {[1e5,1e4],[1e4,1e5]};
% final cooldown
perm_opts.final.iters = 55000;
perm_opts.final.step = 10;
perm_opts.final.temp = [1e6,1e6];

% tuning riff
%{
trials = 5;
[~,~,stat_finals,stats] = corrperm_ampdel_tempering(margs_sort,new_samples,trials,temparams);
corrperm_display_stats(stats);
figure; % save plot of trial run schedule
corrperm_display_stats(stats);
save_current_figure(fullfile(pwd,'ll.tuning_stats'),{'png','pdf'});
mean_err = mean(vertcat(stat_finals{:}))
std_err = std(vertcat(stat_finals{:}))
%}

% command for submitting permutations to LSF
%{
submit_perms('corrperm_ampdel_tempering_module','lsf.submit',ref_dir,perm_dir,100,50,perm_opts);
%}

%% define analysis options
anal_opts = struct;
anal_opts.sig_thresh = 0.95;
anal_opts.power_thresh = 0.1;
anal_opts.tail = 'both';
anal_opts.split_eq = false;
anal_opts.ext = 'll_pair';
anal_opts.lineage_out = unique(H.pcname);
anal_opts.perm_file_mask = 'idx_cell*.mat';
anal_opts.pcount = 1;

% command for running low-level analysis
%{
corrperm_analyze_pairs2(ref_dir,perm_dir,results_dir,anal_opts);
%}
