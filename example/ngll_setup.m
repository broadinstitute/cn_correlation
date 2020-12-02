% script to define global variables for "low level" analysis of 2013 Nature Genetics pan-cancer study

clear

ext = 'ngll'; % distinguishing file extension

ref_dir = fullfile(pwd,'ngll_work/');     % where permutation inputs are stored
perm_dir = fullfile(ref_dir,'ll_permout'); % where permutation output chunk files are stored
results_dir = fullfile(ref_dir,'/ll_results/'); % where significance results are stored


% load disruption profile for permutations
H = import_xrupt('ngll.xrupt.txt');
% load event map for analysis
E = load_emap('ngll.emap');

%% set parameters for low-level tempered annealing
perm_opts = struct;
% misc
perm_opts.algorithm = 'ampdel_tempering';
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

%% define analysis options
anal_opts = struct;
anal_opts.sig_thresh = 0.95;
anal_opts.power_thresh = 0.1;
anal_opts.tail = 'both';
anal_opts.split_eq = false;
anal_opts.ext = ext;
anal_opts.lineage_out = unique(H.pcname);
anal_opts.perm_file_mask = 'idx_cell*.mat';
anal_opts.pcount = 1;
