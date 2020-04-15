% validation code
% recapitulate Travis's correlation analysis for version control factoring 
% (adapted from Nature Genetics submission main correlation script 
% /xchip/gistic/Travis/Correlations/LSF2/0816/NG_ll_wrap.m)
%
% calls corperm_prep to prepare data for 'low level' permutations

clear
dbstop if error
set_verbose_level(40)

% local working directory
reference_dir = [pwd,'/work_dir/'];

% data used in the pancan NG paper
dd = '/xchip/gistic/tcga_projects/pancancer/data_analysis_two/';
gis_dir = [dd,'gistic_analyses/blen_cutoff_series_130814/results/blen0.5_narmp/'];
load([gis_dir,'peak_regs.blen0.5_narmp.mat']); %regs file

% add peak name to regs
amp_del = {'amp','del'};
for k=1:2
    peak_info = read_R_table([dd,'current_dataset/',amp_del{k},'_peak_info.txt']);
    [~,order] = sort([regs{k}.resid_qv]);
    for p=1:length(regs{k})
        regs{k}(order(p)).name = peak_info(p).name;
    end
end

D = load_D([gis_dir,'D.cap1.5.blen0.5_narmp.mat']); %D_file
SIF = '/xchip/gistic/tcga_projects/pancancer/data_analysis_two/current_dataset/archive/pancan12.sample_info_w_subtypes.130815.txt';
sif = read_R_table(SIF);    

sif_names = {sif.tcga_id};

[~,i1,i2] = intersect(sif_names,D.sdesc);

D = reorder_D_cols(D,i2);

%remove "normal"
temp = strmatch('Normal',{sif.subtype});
for i = 1:length(temp)
    sif(temp(i)).subtype = [];
end

%replace "READ" with COAD 
temp = strmatch('READ',{sif.disease});
for i = 1:length(temp)
    sif(temp(i)).disease = 'COAD';
end
D.sis = sif(i1);

% further qualify disease with subtype
for i = 1:length(D.sdesc)
    D.sis(i).disease = [D.sis(i).disease,'_',D.sis(i).subtype];
end

% set the options for low-level analysis
options = struct;
options.event_thresh = 0;           % 0.95 for high-level
options.hilevel = false;            % true for high level
options.minclass_samples = 17;      % 15 for high-level
options.t_amp = 0.3;                % 0.2 for high-level
options.t_del = 0.3;                % 0.2 for high level
options.broad_len_cutoff = 0.50;    % 0.55 for high-level
options.max_disruption = [Inf,Inf]; % [160000,600000] for high level
options.permclass_sisfield = 'disease';

% prepare data for permutations
[D,margs_sort,new_samples] = corrperm_prep(D,regs,reference_dir,options);

%% set parameters for low-level tempered annealing
temparams = struct;
% misc
temparams.algorithm = 'basic_tempered_annealing';
temparams.stat_interval = 1000;
temparams.class_subiters = 100;
temparams.ramp_steps = 10;
temparams.minbinmin = 0;
% initial cycle
temparams.initial.iters = 3e4;
temparams.initial.step = 10;
temparams.initial.temp = [1e5,1e4];
% single tempering cycle
temparams.cycle(1).iters = 5000;
temparams.cycle(1).step = 10;
temparams.cycle(1).heatup = {[1e5,1e4],[1e4,1e5]};
% final cooldown
temparams.final.iters = 55000;
temparams.final.step = 10;
temparams.final.temp = [1e6,1e6];

%% test/tuning runs
trials = 4;
[rand_margs_cell,idx_cell,stat_finals,stats] = corrperm_ampdel_tempering(margs_sort,new_samples,trials,temparams);
corrperm_display_stats(stats);

mean_err = mean(vertcat(stat_finals{:}))
std_err = std(vertcat(stat_finals{:}))

keyboard

%% LSF runs

corrperm_lsf_submission('corrperm_ampdel_tempering_module',ref_dir,perm_dir,500,100,temparams);

%% wait for all those permutations to complete!
keyboard
