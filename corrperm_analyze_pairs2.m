function corrperm_analyze_pairs2(ref_dir,perm_dir,save_dir,options)
% memory-conserving version of corrperm_analyze_pairs

% adapted from /xchip/gistic/Travis/Correlations/LSF2/0816/cancat_mem_saver_and_p.m

%% process optional arguments
if ~exist('options','var')
    options = struct;
end
options = impose_default_value(options,'ext','');           % extension label for output files
options = impose_default_value(options,'sig_thresh',0.25);  % cutoff for q-value
options = impose_default_value(options,'power_thresh',0.1); % cutoff for minimum possible p-value 
options = impose_default_value(options,'split_eq',false);   % set to split observed===permuted between tail and body 
options = impose_default_value(options,'pcount',1);         % pseudo count (1 or 0) 
options = impose_default_value(options,'perm_file_mask','idx_cell*');% set to split obs==permuted between tail and body 

if options.split_eq
     warning('The ''split_eq'' way of calculating P-values is deprecated.');
     options.split_eq = false;
end

% fix output extension
ext = options.ext;
if ~isempty(ext) && ext(1) ~= '.'
    ext = ['.',ext];
end

% create results directory, if necessary
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

%% read reference inputs
verbose('Loading inputs from reference directory ''%s''',10,ref_dir);
[~,Binary,Binary_amps,Binary_dels,samples,regs] = load_permutation_ref_inputs(ref_dir);
% Binary are the observed events IN ORIGINAL UNOPTIMIZED SAMPLE ORDER
%   (logical matrix, Nevents x Nsamples)
% Binary_amps, Binary_dels is Binary separated into amps and dels
% samples is the lineage partition, with index order optimized
% regs are peak regions (with 'name' field added) 

% lineage output controlled by a list of strings aligned with 'samples'
% partition into lineages - if empty, no output
if ~isfield(options,'lineage_out') || isempty(options.lineage_out)
    % default lineage output is none
    options.lineage_out = cell(size(samples));
end

Nlineages = length(samples);

Namps = length(regs{1});
Ndels = length(regs{2}); 

chr_amps = [regs{1}.chrn];
chr_dels = [regs{2}.chrn];

[Nevents,Nsamples] = size(Binary);

%% calculate observed co-occurrences
Binary = Binary+0;
% by lineage
lin_obs = zeros(Nevents,Nevents,Nlineages);
for i = 1:Nlineages
    lin_obs(:,:,i) = Binary(:,samples{i})*Binary(:,samples{i})';
end

% overall
obs_tot = sum(lin_obs,3);


%% storage for statisitics
% (in this version, we calculate these as we pass over the files)

% overall
le_counts = zeros(Nevents,Nevents);
ge_counts = zeros(Nevents,Nevents);
eq_counts = zeros(Nevents,Nevents);

% by lineage
le_counts_byclass = zeros(Nevents,Nevents,Nlineages);
ge_counts_byclass = zeros(Nevents,Nevents,Nlineages);
eq_counts_byclass = zeros(Nevents,Nevents,Nlineages);

Nperms = 0; % initialize permutation count

%% loop over permutation results files

files = dir([perm_dir,options.perm_file_mask]);

verbose('Reading %d permutation chunks from ''%s''',10,length(files),perm_dir);

for k = 1:length(files)
    verbose(files(k).name,10);
    tic
    load([perm_dir,files(k).name]);  % 'idx_cell' cell array, each element NCHR x Nsamples
    npf = length(idx_cell);

    % get number of chromosomes from first chunck of permutations
    if ~exist('NCHR','var')
        NCHR = size(idx_cell{1},1);
        % initialize chromosome-structure storage
        chrns_a = cell(1,NCHR);
        chrns_d = cell(1,NCHR);
        for i = 1:NCHR
            chrns_a{i} = find(chr_amps==i);
            chrns_d{i} = find(chr_dels==i);
        end

    end

    % loop over permutations in the chunk
    for i = 1:npf
        Nperms = Nperms+1;
        % map chromosome-swap indices to event swaps
        Bi_a = false(Namps,Nsamples);
        Bi_d = false(Ndels,Nsamples);
        % loop over chromosomes in the permutation 
        for j = 1:NCHR
            % use permutation indices to rearrange observed events
            idx_mat = idx_cell{i};
            % check for old 3D index matrix and convert if necessary
            if length(size(idx_mat))==3
                idx_mat = idx_mat(:,:,1);
            end
            Bi_a(chrns_a{j},:) = Binary_amps(chrns_a{j},idx_mat(j,:));
            Bi_d(chrns_d{j},:) = Binary_dels(chrns_d{j},idx_mat(j,:));
        end
        % map events and their co-occurrances
        emap = double([Bi_a;Bi_d]); % event map for permutation (Nevents X Nsamples)
        cooc = emap*emap'; % co-ocurrence count for permutation (Nevents X Nevents)
        le_counts = le_counts + (cooc <= obs_tot);
        ge_counts = ge_counts + (cooc >= obs_tot);
        eq_counts = eq_counts + (cooc == obs_tot);
        % also accumulate per-lineage le & ge counts
        for l = 1:Nlineages
            submap = emap(:,samples{l});
            lcooc = submap*submap';
            le_counts_byclass(:,:,l) = le_counts_byclass(:,:,l) + (lcooc <= lin_obs(:,:,l));
            ge_counts_byclass(:,:,l) = ge_counts_byclass(:,:,l) + (lcooc >= lin_obs(:,:,l));
            eq_counts_byclass(:,:,l) = eq_counts_byclass(:,:,l) + (lcooc == lin_obs(:,:,l));
        end
    end
    toc
end

%% count events for lineage-specfic correlation and anti-correlation p-values

% event => chromosome map
chrns = [chr_amps';chr_dels'];

% calculate number of pairs on different chromosomes
Npairs = 0;
for i = 1:length(chrns)
    Npairs = Npairs + sum(chrns(i) < chrns);
end

pscnt = options.pcount;

% allocate storage for lineage-specific p-values
verbose('calculating p-values for lineage specific co-occurrences',10);
tic
regs_idx = zeros(Npairs,2);            % list of unique event pair indices
p_list_corr = zeros(Npairs,Nlineages); % correlation p-values per lineage
p_list_anti = zeros(Npairs,Nlineages); % anti-correlation p-values per lineage
for l = 1:Nlineages
    s = 1; % pair index
    p_list_cpow = zeros(Npairs,1);  % min fisher exact p-value for correlation
    p_list_apow = zeros(Npairs,1);  % min fisher exact p-value for anti-correlation
   for i = 1:Nevents
        for j = i+1:Nevents
            if chrns(i)~=chrns(j)
                % maximum power calculation
                [p_list_cpow(s),p_list_apow(s)] = max_fish_power(length(samples{l}), sum(Binary(i,samples{l})), sum(Binary(j,samples{l}))); 
                % two ways of dealing with permuted==observed
                if options.split_eq
                    p_list_corr(s,l) = (ge_counts_byclass(i,j,l) - eq_counts_byclass(i,j,l)/2 + pscnt) / (Nperms + pscnt); 
                    p_list_anti(s,l) = (le_counts_byclass(i,j,l) - eq_counts_byclass(i,j,l)/2 + pscnt) / (Nperms + pscnt);
                else
                    p_list_corr(s,l) = (ge_counts_byclass(i,j,l) + pscnt) / (Nperms + pscnt);
                    p_list_anti(s,l) = (le_counts_byclass(i,j,l) + pscnt) / (Nperms + pscnt);
                    % (added +1 pseudocount 2016-01-20)
                end
                regs_idx(s,:) = [i,j];
                s = s+1;
            end
        end
    end
    if l <= length(options.lineage_out) && ~isempty(options.lineage_out{l})
        save_pair_pvalues(regs,[p_list_anti(:,l),regs_idx,p_list_apow],...
            [save_dir,'anticorr_pair.',options.lineage_out{l},ext,'.txt'],1,options.sig_thresh,options.power_thresh);
        save_pair_pvalues(regs,[p_list_corr(:,l),regs_idx,p_list_cpow],...
            [save_dir,'correlate_pair.',options.lineage_out{l},ext,'.txt'],1,options.sig_thresh,options.power_thresh);
    end
end
toc

%% count events for overall correlation and anti-correlation p-values

% sum permuted per-lineage event coincidences into overall permuted event coincidences
verbose('counting overall event co-occurrences',10)

% allocate storage for lineage-specific p-values
regs_idx = zeros(Npairs,2);   % peak indices
p_corr = zeros(Npairs,1);
p_anti = zeros(Npairs,1);
p_cpow = zeros(Npairs,1);  % min fisher exact p-value for correlation
p_apow = zeros(Npairs,1);  % min fisher exact p-value for anti-correlation
% calculate p-values
verbose('calculating overall p-values for co-occurrences',10);
tic
s = 1; % pair index
for i = 1:Nevents-1
    for j = i+1:Nevents
        if chrns(i)~=chrns(j)
            %% maximum power calculation
            [p_cpow(s),p_apow(s)] = max_fish_power(Nsamples,sum(Binary(i,:)),sum(Binary(j,:)));
            %% p value calculations
            if options.split_eq
                p_corr(s) = (ge_counts(i,j) - eq_counts(i,j)/2 + pscnt) / (Nperms + pscnt);
                p_anti(s) = (le_counts(i,j) - eq_counts(i,j)/2 + pscnt) / (Nperms + pscnt); 
            else
                p_corr(s) = (ge_counts(i,j) + pscnt) / (Nperms + pscnt);
                p_anti(s) = (le_counts(i,j) + pscnt) / (Nperms + pscnt);
            end
            % save pair index 
            regs_idx(s,:) = [i,j];
            s = s+1;
        end
    end
end
toc

%% save results files
verbose('saving results',10);toc

% save binaries for forensics
save([save_dir 'pair_results',ext,'.mat'],'regs_idx','p_corr','p_anti','p_cpow','p_apow');
%!save([save_dir 'pair_perm_tot',ext,'.mat'],'perm_tot');
save([save_dir 'pair_obs_tot',ext,'.mat'],'obs_tot');

% save text results files
save_pair_pvalues(regs,[p_anti,regs_idx,p_apow],[save_dir,'anticorr_pair',ext,'.txt'],...
                1,options.sig_thresh,options.power_thresh);
save_pair_pvalues(regs,[p_corr,regs_idx,p_cpow],[save_dir,'correlate_pair',ext,'.txt'],...
                1,options.sig_thresh,options.power_thresh);

