% script to define global variables for "low level" analysis of 2013 Nature Genetics pan-cancer study

clear
warning('UNVALIDATED NGHL EXAMPLE'); %!!! validate some day...

ext = 'nghl'; % distinguishing file extension
ref_dir = fullfile(pwd,'hl_work/'); 
perm_dir = fullfile(ref_dir,'hl_permout'); % where permutation outputs are stored
results_dir = fullfile(ref_dir,'/hl_results/');

% load disruption profile for permutations
H = import_xrupt('nghl.xrupt.txt');
% load event map for analysis
E = load_emap('nghl.emap');

%% set parameters for high-level tempered annealing
perm_opts = struct;

% misc
perm_opts.algorithm = 'ampdel_tempering';
perm_opts.stat_interval = 1000;
perm_opts.class_subiters = 100;
perm_opts.ramp_steps = 10;
perm_opts.minbinmin = 100;
% initial cycle
perm_opts.initial.iters = 8e4;
perm_opts.initial.step = 10;
perm_opts.initial.temp = [50,1e3];
% 5 tempering cycles
for i = 1:5
    perm_opts.cycle(i).iters = 15000;
    perm_opts.cycle(i).step = 10;
    perm_opts.cycle(i).heatup = {[1e6,5e3]*(10^i),[5e3,1e6]*(10^i)};
end
% final cooldown
perm_opts.final.iters = 55000;
perm_opts.final.step = 10;
perm_opts.final.temp = [1e6,1e6];

% define analysis options
anal_opts = struct;
anal_opts.sig_thresh = 0.95;
anal_opts.power_thresh = 0.1;
anal_opts.tail = 'both';
anal_opts.split_eq = false;
anal_opts.ext = ext;
anal_opts.lineage_out = unique({D.sis.disease});
anal_opts.perm_file_mask = 'idx_cell*.mat';
anal_opts.pcount = 1;
