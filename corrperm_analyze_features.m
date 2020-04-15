function corrperm_analyze_features(features,ref_dir,perm_dir,save_dir,options)
%CORRPERM_ANALYZE_FEATURES analyse permutation results for crrelatio with sample features
%   corrperm_analyze_features(FEATURES,REF_DIR,PERM_DIR,SAVE_DIR,OPTIONS)
%
% FEATURES - a list of supplementary data fields in D (matching D.supacc)
%           to use as sample features
% REF_DIR - a directory containing several data files:
%               Binary_amps.mat - matrix of of observed amp.events x samples
%               Binary_dels.mat - matrix of of observed del.events x samples
%               new_samples.mat - defines partition of samples into lineages
%               peak_regs.mat - a GISTIC regs struct to wch has been added
%                       a 'name' field for the output
%               D.mat - D-struct (required, but not used by this function)
% PERM_DIR - directory where permutation results were stored
% SAVE_DIR - directory where results from this function will be saved
% OPTIONS - optional parameters:
%   OPTIONS.ext - results file extenion (default none)
%!!! finsh documenting
%

% sketch: function to correlate peaks with sample classes 
%!!! what will be function inputs
%{
ref_dir = '/xchip/gistic/schum/mcmc_correlation_test_130930/work_dir/'; % prepared inputs
perm_dir = '/xchip/gistic/schum/mcmc_correlation_test_130930/output_dir/'; % permutations
save_dir = '/xchip/gistic/schum/mcmc_correlation_test_130930/results_wgd/';  %  final results (output)

ext = 'll_0.3_sub'; % extension for output files
%}

% adapted from /xchip/gistic/Travis/Correlations/LSF2/0816/cancat_mem_saver_and_p.m

if ~exist('options','var')
    options = struct;
end
options = impose_default_value(options,'ext','');
options = impose_default_value(options,'sig_thresh',0.25);
options = impose_default_value(options,'power_thresh',0.1); % cutoff for minimum possible p-value 
options = impose_default_value(options,'two_tailed',false);
options = impose_default_value(options,'tail','both',{'both','corr','anti'});
options = impose_default_value(options,'fit_quantile',1.00);

% fix output extension
ext = options.ext;
if ~isempty(ext) && ext(1) ~= '.'
    ext = ['.',ext];
end
options.ext = ext;

% read reference inputs
verbose('Loading inputs from reference directory ''%s''',10,ref_dir);
[D,Binary,Binary_amps,Binary_dels,samples,regs] = load_permutation_ref_inputs(ref_dir);
% Binary are the observed events (logical matrix, Nevents x Nsamples)
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
[Nevents,Nsamples] = size(Binary);

%% create binary matrix of features
feature_matrix = D.supdat ~= 0;

if ~isempty(features)
    if ischar(features)
        features = {features};
    end
    keeper_features = ismember(cellstr(D.supacc),features);
    feature_matrix = feature_matrix(keeper_features,:);
end
Nfeatures = size(feature_matrix,1);

Namps = length(regs{1});
Ndels = length(regs{2}); 


% chromosome locations of peaks
chr_amps = [regs{1}.chrn];
chr_dels = [regs{2}.chrn];


%% Actually grab output and analyze permutations
files = dir([perm_dir,'idx_cell*']);

Nperms = 0; % overall permutation index

verbose('Reading permutations from ''%s''',10,perm_dir);
% loop over permutation results files
for k = 1:length(files)
    verbose(files(k).name,10);
    load([perm_dir,files(k).name]);  % 'idx_cell', each element NCHR x Nsamples x amp/del
    npf = length(idx_cell);
    % create storage for permutation results across files
    if ~exist('Binary1','var')
        Binary1 = false(Nevents,Nsamples,npf*length(files));
    end
    % get number of chromosomes from permutations
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
    Bi_a = cell(1,npf);
    Bi_d = cell(1,npf);
    % loop over permutations in the file
    for i = 1:npf
        Bi_a{i} = false(Namps,Nsamples);
        Bi_d{i} = false(Ndels,Nsamples);
        % loop over chromosomes in the permutation
        for j = 1:NCHR
            % use permutation indices to rearrange observed events
            idx_mat = idx_cell{i};
            % check for old 3D index matrix and convert if necessary
            if length(size(idx_mat))==3
                idx_mat = idx_mat(:,:,1);
            end
            Bi_a{i}(chrns_a{j},:) = Binary_amps(chrns_a{j},idx_mat(j,:));
            Bi_d{i}(chrns_d{j},:) = Binary_dels(chrns_d{j},idx_mat(j,:));
        end
    end
    %!!! lin_couns not used
    %{
    % initialize lineage counts
    if ~exist('lin_couns','var')
        lin_couns = cell(1,Nlineages);
        for i = 1:Nlineages
            lin_couns{i} = zeros(Nevents,Nfeatures,npf*length(files));
        end
    end
    %}

    % count co-occurrances
    for i = 1:npf
        Nperms = Nperms+1;
        % accumulate events in Nevents X Nsamples X Nperms matrix
        Binary1(:,:,Nperms) = [Bi_a{i};Bi_d{i}]; 
        %!!! lin_couns not used
        %{
        % also accumulate per-lineage
        for l = 1:Nlineages
            lin_couns{l}(:,:,Nperms) = double(Binary1(:,samples{l},Nperms)) * double(feature_matrix(:,samples{l}))';
        end
        %}
    end
end

%% if we are being selective about the fits, reduce permutations to those with the best fit 
if options.fit_quantile < 1.00
    [afit,dfit] = corrperm_get_final_stats([perm_dir,'stat_finals.*.mat']);
    rms_fit = sqrt(afit.^2+dfit.^2);
    [~,order] = sort(rms_fit);
    keepers = order <= round(options.fit_quantile*length(order));
    Binary1 = Binary1(:,:,keepers);
    Nperms = size(Binary1,3);
    verbose('keeping %d best-fitting permutations',20,Nperms);
end
      
% tabulate observed event incidence with feature per lineage
Binary = Binary+0;
lin_obs = cell(1,Nlineages);
for i = 1:Nlineages
    lin_obs{i} = Binary(:,samples{i}) * feature_matrix(:,samples{i})';
end

obs_tot = zeros(size(lin_obs{1}));

% sum lineage-specific event coincidences into overall event coincidences
for i = 1:Nlineages
    obs_tot = obs_tot + lin_obs{i};
end

perm_tot = zeros(Nevents,Nfeatures,Nperms);
for i = 1:Nperms
    perm_tot(:,:,i) = double(Binary1(:,:,i)) * double(feature_matrix');
end


%% count events vs features for correlation and anti-correlation p-values

%% save results files
verbose('saving results',10)

% create results directory, if necessary
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end
%!R = p_power_count(Binary,perm_tot,obs_tot,feature_matrix);
R = analyze_and_save(save_dir,Binary,perm_tot,obs_tot,feature_matrix,regs,features,options,ext);

%{
% calculate number of pairs of events versus features
Npairs = Nevents*Nfeatures;

% allocate storage for lineage-specific p-values
regs_idx = zeros(Npairs,2);   % peak indices
p_corr = zeros(Npairs,1);
p_anti = zeros(Npairs,1);
% variables for observational summary an,regsd power
p_cpow = zeros(Npairs,1);  % min fisher exact p-value for correlation
p_apow = zeros(Npairs,1);  % min fisher exact p-value for anti-correlation

% calculate p-values
s = 1; % pair index
for i = 1:Nevents
    for j = 1:Nfeatures
        % maximum power calculation
        [p_cpow(s),p_apow(s)] = max_fish_power(Nsamples,sum(Binary(i,:)),sum(feature_matrix(j,:)));

        % p value calculation
        p_corr(s) = (sum(squeeze(perm_tot(i,j,:))>=obs_tot(i,j)))./Nperms;
        p_anti(s) = (sum(squeeze(perm_tot(i,j,:))<=obs_tot(i,j)))./Nperms;
        regs_idx(s,:) = [i,j];
        s = s+1;
    end
end
%}

% save binaries for forensics
save([save_dir 'feature_results',ext,'.mat'],'R');
save([save_dir 'feature_perm_tot',ext,'.mat'],'perm_tot','-v7.3');
save([save_dir 'feature_obs_tot',ext,'.mat'],'obs_tot');

%{
% save text depending on user's choice of tail
if strcmp(options.tail,'both') || options.two_tailed
    % correct for two-tailed test
    bothlist = [ R.p_corr,R.regs_idx,R.p_cpow,ones(size(R.p_corr));
                R.p_anti,R.regs_idx,R.p_apow,zeros(size(R.p_anti)) ];
    save_feature_pvalues(regs,features,bothlist,...
            [save_dir,'peak_vs_feature',ext,'.txt'],1,options.sig_thresh);
        
elseif strcmp(options.tail,'corr')
    % save just correlations
    save_feature_pvalues(regs,features,[R.p_corr,R.regs_idx,R.p_cpow],...
                [save_dir,'correlate_peak_vs_feature',ext,'.txt'],1,options.sig_thresh);
    % do individual tails (!!!NOTE should chose one with an option)
elseif strcmp(options.tail,'anti')
    % save just anticorrelations
    save_feature_pvalues(regs,features,[R.p_anti,R.regs_idx,R.p_apow],...
                [save_dir,'anticorr_peak_vs_feature',ext,'.txt'],1,options.sig_thresh);
else
    % invalid optionCHR
    throw(MException('snp_correlation:bad_option',...
                     '''%s'' is not a valid parameter',options.tail));
end
%}

%% lineage-specific correlations

linouts = find(~cellfun(@isempty,options.lineage_out));
for l = 1:length(linouts)
    lx = linouts(l);
    lintext = options.lineage_out{lx};
    lperm_tot = zeros(Nevents,Nfeatures,Nperms);
    
    % sum up the permutations for this lineage
    for i = 1:Nperms
        lperm_tot(:,:,i) = double(Binary1(:,samples{lx},i)) * double(feature_matrix(:,samples{lx})');
    end
    analyze_and_save(save_dir,Binary(:,samples{lx}),lperm_tot,lin_obs{lx},feature_matrix(:,samples{lx}),...
                    regs,features,options,['.',lintext,ext]);
%{
    % calculate power and maximum p-value
    R = p_power_count(Binary(:,samples{lx}),lperm_tot,lin_obs{lx},feature_matrix(:,samples{lx}));
    % save text depending on user's choice of tail
    if strcmp(options.tail,'both') || options.two_tailed
        % correct for two-tailed test
        bothlist = [ R.p_corr,R.regs_idx,R.p_cpow,ones(size(R.p_corr));
                    R.p_anti,R.regs_idx,R.p_apow,zeros(size(R.p_anti)) ];
        save_feature_pvalues(regs,features,bothlist,...
                [save_dir,'peak_vs_feature.',lintext,ext,'.txt'],1,options.sig_thresh);

    elseif strcmp(options.tail,'corr')
        % save just correlations
        save_feature_pvalues(regs,features,[R.p_corr,R.regs_idx,R.p_cpow],...
                    [save_dir,'correlate_peak_vs_feature.',lintext,ext,'.txt'],1,options.sig_thresh);
        % do individual tails (!!!NOTE should chose one with an option)
    elseif strcmp(options.tail,'anti')
        % save just anticorrelations
        save_feature_pvalues(regs,features,[R.p_anti,R.regs_idx,R.p_apow],...
                    [save_dir,'anticorr_peak_vs_feature.',lintext,ext,'.txt'],1,options.sig_thresh);
    else
        % invalid option
        throw(MException('snp_correlation:bad_option',...
                         '''%s'' is not a valid parameter',options.tail));
    end
%}
end

%% subfunction: save results
function R = analyze_and_save(save_dir,Binary,perm_tot,obs_tot,feature_matrix,regs,features,options,ext)
    R = p_power_count(Binary,perm_tot,obs_tot,feature_matrix);
    % save text depending on user's choice of tail
    if strcmp(options.tail,'both') || options.two_tailed
        % correct for two-tailed test
        bothlist = [ R.p_corr,R.regs_idx,R.p_cpow,ones(size(R.p_corr));
                    R.p_anti,R.regs_idx,R.p_apow,zeros(size(R.p_anti)) ];
        save_feature_pvalues(regs,features,bothlist,...
                [save_dir,'peak_vs_feature',ext,'.txt'],1,options.sig_thresh);

    elseif strcmp(options.tail,'corr')
        % save just correlations
        save_feature_pvalues(regs,features,[R.p_corr,R.regs_idx,R.p_cpow],...
                    [save_dir,'correlate_peak_vs_feature',ext,'.txt'],1,options.sig_thresh);
        % do individual tails (!!!NOTE should chose one with an option)
    elseif strcmp(options.tail,'anti')
        % save just anticorrelations
        save_feature_pvalues(regs,features,[R.p_anti,R.regs_idx,R.p_apow],...
                    [save_dir,'anticorr_peak_vs_feature',ext,'.txt'],1,options.sig_thresh);
    else
        % invalid option
        throw(MException('snp_correlation:bad_option',...
                         '''%s'' is not a valid parameter',options.tail));
    end


%% subfunction: count events vs features for correlation and anti-correlation p-values
function R = p_power_count(Binary,perm_tot,obs_tot,feature_matrix)

[Nevents,Nsamples,~] = size(Binary);
Nfeatures = size(feature_matrix,1);
Nperms = size(perm_tot,3);

% calculate number of pairs of events versus features
Npairs = Nevents*Nfeatures;

R = struct;
% allocate storage for lineage-specific p-values
R.regs_idx = zeros(Npairs,2);   % peak indices
R.p_corr = zeros(Npairs,1);
R.p_anti = zeros(Npairs,1);
% variables for observational summary and power
R.p_cpow = zeros(Npairs,1);  % min fisher exact p-value for correlation
R.p_apow = zeros(Npairs,1);  % min fisher exact p-value for anti-correlation

% fisher exact p-values
R.p_fisher = zeros(Npairs,2); % actual FET results (each row [left-tailed,right-tailed])
% event counts
R.peak_count = zeros(Npairs,1);
R.feat_count = zeros(Npairs,1);
R.co_occurrances = zeros(Npairs,1);

% calculate p-values
s = 1; % pair index
for i = 1:Nevents
    for j = 1:Nfeatures
        % maximum power calculation
        peakmarg = sum(Binary(i,:));
        featmarg = sum(feature_matrix(j,:));
        [R.p_cpow(s),R.p_apow(s)] = max_fish_power(Nsamples,peakmarg,featmarg);
        
        ad = sum(Binary(i,:)==feature_matrix(j,:));
        a = (Nsamples+ad-peakmarg-featmarg)/2;
        d = ad - a;
        b = featmarg - d;
        c = peakmarg - d;
        R.p_fisher(s,1) = fisher_exact_test(a,b,c,d);
        R.p_fisher(s,2) = fisher_exact_test(a,b,c,d,[],'left');
%!        [R.p_fisher(s,1),R.p_fisher(s,2)] = fisher_exact_test(a,b,c,d);
        R.peak_count(i) = peakmarg;
        R.feat_count(i) = featmarg;
        R.co_occurrances(i) = d;

        % p value calculation
        R.p_corr(s) = (sum(squeeze(perm_tot(i,j,:))>=obs_tot(i,j)))./Nperms;
        R.p_anti(s) = (sum(squeeze(perm_tot(i,j,:))<=obs_tot(i,j)))./Nperms;
        R.regs_idx(s,:) = [i,j];
        s = s+1;
    end
end


