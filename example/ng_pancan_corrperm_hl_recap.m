% calls corrperm_prep to prepare data for 'high level' permutations

clear
%dbstop if error
set_verbose_level(40)

% local working directories
%ref_dir = [pwd,'/hl_work/'];     % where permutation inputs are stored
ref_dir = '/xchip/beroukhimlab/gistic/corrperm/hl_work/';
perm_dir = fullfile(ref_dir,'/hl_permout/'); % where permutation outputs are stored
results_dir = fullfile(ref_dir,'/hl_results/');

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

% set the options for preparing high-level analysis
options = struct;
options.hilevel = true;             % false for low-level
options.event_thresh = 0.95;        % 0 for high-level
options.minclass_samples = 17;      % 17 for low-level
options.t_amp = 0.2;                % 0.3 for low-level
options.t_del = 0.2;                % 0.3 for low-level
options.broad_len_cutoff = 0.55;    % 0.50 for low-level
options.max_disruption = [160000,600000];  % [Inf,Inf] for low-level 
options.permclass_sisfield = 'disease';


% prepare data for permutations
[D,margs_sort,new_samples] = corrperm_prep(D,regs,ref_dir,options);

%% set parameters for low-level tempered annealing
temparams = struct;
% misc
temparams.algorithm = 'basic_tempered_annealing';
temparams.stat_interval = 1000;
temparams.class_subiters = 100;
temparams.ramp_steps = 10;
temparams.minbinmin = 100;
% initial cycle
temparams.initial.iters = 8e4;
temparams.initial.step = 10;
temparams.initial.temp = [50,1e3];
% 5 tempering cycles
for i = 1:5
    temparams.cycle(i).iters = 15000;
    temparams.cycle(i).step = 10;
    temparams.cycle(i).heatup = {[1e6,5e3]*(10^i),[5e3,1e6]*(10^i)};
end
% final cooldown
temparams.final.iters = 55000;
temparams.final.step = 10;
temparams.final.temp = [1e6,1e6];

%% test/tuning runs
trials = 4;
[~,~,stat_finals,stats] = corrperm_ampdel_tempering(margs_sort,new_samples,trials,temparams);
keyboard

% save plot of trial run schedule
figure;
corrperm_display_stats(stats);
save_current_figure(fullfile(pwd,'hl.tuning_stats'),{'png','pdf'});

mean_err = mean(vertcat(stat_finals{:}))
std_err = std(vertcat(stat_finals{:}))

%! no tuning pause for now
%! keyboard

%% LSF runs

submit_perms('corrperm_ampdel_tempering_module','lsf.submit',ref_dir,perm_dir,100,50,temparams);

% wait for all those permutations to complete!
keyboard

% analyze pair correlations
options = struct;
options.sig_thresh = 0.95;
options.power_thresh = 0.1;
options.tail = 'both';
options.split_eq = true;
options.ext = 'hl_pair';
options.lineage_out = unique({D.sis.disease});
options.perm_file_mask = 'idx_cell.chunk.*.mat';
options.pcount = 0;

corrperm_analyze_pairs2(ref_dir,perm_dir,results_dir,options);
