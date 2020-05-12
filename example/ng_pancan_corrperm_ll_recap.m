% validation code
% recapitulate Travis's correlation analysis for version control factoring
% (adapted from Nature Genetics submission main correlation script 
% /xchip/gistic/Travis/Correlations/LSF2/0816/NG_ll_wrap.m)
%
% calls corperm_prep to prepare data for 'low level' permutations

clear
%dbstop if error
set_verbose_level(40)

% local working directories
ref_dir = fullfile(pwd,'ll_work/');     % where permutation inputs are stored
perm_dir = fullfile(ref_dir,'ll_permout'); % where permutation output chunk files are stored
results_dir = fullfile(ref_dir,'/ll_results/'); % where significance results are stored
% input files
gistic_peaks_file = fullfile(pwd,'peak_regs.blen0.5_narmp.mat'); % peak results from GISTIC
gistic_data_file = fullfile(pwd,'/D.cap1.5.blen0.5_narmp.mat');   % copy-number data w/ziggurat 
peak_name_files = {fullfile(pwd,'/amp_peak_info.txt'),fullfile(pwd,'/del_peak_info.txt')}; % peak friendly names
subtype_info_file = fullfile(pwd,'/pancan12.sample_info_w_subtypes.130815.txt');   % subtype information

% "bare load" matlab peak file into variable 'regs'
load(gistic_peaks_file);

% add manually edited peak name to regs (for results)
amp_del = {'amp','del'};
for k=1:2
    peak_info = read_R_table(peak_name_files{k});
    [~,order] = sort([regs{k}.resid_qv]);
    for p=1:length(regs{k})
        regs{k}(order(p)).name = peak_info(p).name;
    end
end


%% Fix sample information in D.sis
% this analysis used subtypes to further qualify lineage / restrict permutations

%% type/subtype processing
sif = read_R_table(subtype_info_file);    

% blank "Normal" subtypes
temp = strmatch('Normal',{sif.subtype});
for i = 1:length(temp)
    sif(temp(i)).subtype = '';
end

% treat rectal as colon adenocarcinoma: replace "READ" with COAD 
temp = strmatch('READ',{sif.disease});
for i = 1:length(temp)
    sif(temp(i)).disease = 'COAD';
end

% load copy number data and update sample information
D = load_D(gistic_data_file);
[~,i1,i2] = intersect({sif.tcga_id},D.sdesc);
D = reorder_D_cols(D,i2);
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

%% prepare data for permutations
[D,margs_sort,new_samples] = corrperm_prep(D,regs,ref_dir,options);

%% set parameters for low-level tempered annealing
temparams = struct;
% misc
temparams.algorithm = 'basic_tempered_annealing'; % (for validation, does not set algorithm)
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

% save trial run schedule
figure;
corrperm_display_stats(stats);
save_current_figure(fullfile(pwd,'ll.tuning_stats'),{'png','pdf'});

mean_err = mean(vertcat(stat_finals{:}))
std_err = std(vertcat(stat_finals{:}))

%tuning pause
keyboard

%% LSF runs
% do 5000 permutations
corrperm_lsf_submission('corrperm_ampdel_tempering_module',ref_dir,perm_dir,100,50,temparams);
%!corrperm_uger_submission('corrperm_ampdel_tempering_module',ref_dir,perm_dir,100,50,temparams);

% wait for all those permutations to complete!
keyboard

%% analyze pair correlations
options = struct;
options.sig_thresh = 0.95;
options.power_thresh = 0.1;
options.tail = 'both';
options.split_eq = true;
options.ext = 'll_pair';
options.lineage_out = unique({D.sis.disease});
options.perm_file_mask = 'idx_cell.chunk.*.mat';
options.pcount = 0;

corrperm_analyze_pairs2(ref_dir,perm_dir,results_dir,options);
%corrperm_analyze_features(features,ref_dir,perm_dir,ref_dir,options);
