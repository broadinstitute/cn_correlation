% script to create "low level" disruption profile and event map files
% from 2013 GISTIC results used in Pan-Cancer 2013 Nat Gen paper

clear
set_verbose_level(40)

ref_dir = fullfile(pwd,'ngll_work/');     % where permutation inputs are stored

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

% prepare data for permutations
[D,margs_sort,pcindex] = corrperm_prep(D,regs,ref_dir,options);

%% save prepared data

% create disruption profile, export to file
H = create_xrupt(D,margs_sort,pcindex);
export_xrupt(H,fullfile(pwd,'ngll.xrupt.txt'));
save(fullfile(pwd,'ngll.xrupt.mat'),'H');

% create event map from corrperm_prep output files, export to files
load(fullfile(ref_dir,'Binary_amps.mat'));
load(fullfile(ref_dir,'Binary_dels.mat'));
E = create_emap(D,regs,pcindex,Binary_amps,Binary_dels);
%!!! format in development
%!!!export_emap(E,fullfile(pwd,'ngll.emap.txt'));
save(fullfile(pwd,'ngll.emap.mat'),'E');
