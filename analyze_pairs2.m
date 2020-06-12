function analyze_pairs2(E,perm_dir,save_dir,options)
% memory-conserving version of corrperm_analyze_pairs
%!!!DOCME

%% process optional arguments
if ~exist('options','var')
    options = struct;
end

% optional argument defaults
options = impose_default_value(options,'ext','');           % extension label for output files
options = impose_default_value(options,'sig_thresh',0.25);  % cutoff for q-value
options = impose_default_value(options,'power_thresh',0.1); % cutoff for minimum possible p-value 
options = impose_default_value(options,'split_eq',false);   % set to split observed===permuted between tail and body 
options = impose_default_value(options,'pcount',1);         % pseudo count (1 or 0) 
options = impose_default_value(options,'perm_file_mask','idx_cell*');
options = impose_default_value(options,'lineage_out',{});
options = impose_default_value(options,'analyze_lineages',false);%
options = impose_default_value(options,'lineage_out',E.pcname);    % filter list of lineage subsets to analyze

options % display

% input warnings
%!if options.split_eq
%!     warning('The ''split_eq'' way of calculating P-values is deprecated.');
%!end
if options.pcount == 0
    warning('A pseudocount of 0 is deprecated.');
end
if options.analyze_lineages && ~isempty(options.lineage_out)
    if isempty(E.pcx)
        warning('no lineage info in event map, canceling lineage analysis');
        options.analyze_lineages = false;
    else
        unklin = ~ismember(options.lineage_out,E.pcname);
        if any(unklin)
            warning(strjoin({'unknown lineages requested:',options.lineage_out{unklin}}));
            if all(unklin)
                warning('canceling lineage analysis');
                options.analyze_lineages = false;
            end
        end
    end
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

% process lineage filter
lindices = {};
linames = {};
if options.analyze_lineages
    % retain matchng lineages
    linx = ismember(E.pcname,options.lineage_out);
    if sum(linx) > 0
        lindices = E.pcx(linx);
        linames = E.pcname(linx);
    end
end
Nlineages = length(lindices);

% get sizes from input data
[Nevents,Nsamples] = size(E.dat);
Nchr = length(unique(E.event.chrn));

%% calculate observed co-occurrences
% overall
emap = double(E.dat);
obs_tot = emap * emap';
% by lineage
if Nlineages > 0
    lin_obs = zeros(Nevents,Nevents,Nlineages);
    for l = 1:Nlineages
        emap = double(E.dat(:,lindices{l}));
        lin_obs(:,:,l) = emap * emap';
    end
end

%% storage for accumulated statisitics
%
% This version of the analysis aggregates three counts for each event pair:
% equality with, fewer than, and greater than the observed co-occurrance.
%
% overall
le_counts = zeros(Nevents,Nevents);
ge_counts = zeros(Nevents,Nevents);
eq_counts = zeros(Nevents,Nevents);
% by lineage
if Nlineages > 0
    le_counts_byclass = zeros(Nevents,Nevents,Nlineages);
    ge_counts_byclass = zeros(Nevents,Nevents,Nlineages);
    eq_counts_byclass = zeros(Nevents,Nevents,Nlineages);
end
%%%%
% chromosome => event mapping
chrns_e = cell(1,Nchr);
for c = 1:Nchr
    chrns_e{c} = find(E.event.chrn == c);
end

Nperms = 0; % initialize permutation count

%% chunk processing

files = dir(fullfile(perm_dir,options.perm_file_mask));
Nfiles = length(files);
verbose('Reading %d permutation chunks from ''%s''',10,Nfiles,perm_dir);
if ~Nfiles
    error('No permutation data to process');
end
for k = 1:Nfiles
    verbose(files(k).name,10);
    tic
    load(fullfile(perm_dir,files(k).name));  % 'idx_cell' cell array, each element Nchr x Nsamples
    npf = length(idx_cell);

    % loop over permutations in the chunk
    for i = 1:npf
        % count permutation
        Nperms = Nperms+1;        

        % use permutation indices to rearrange observed events
        idx_mat = idx_cell{i};
        % check for old 3D index matrix and convert if necessary
        if length(size(idx_mat))==3
            idx_mat = idx_mat(:,:,1);
        end
        % map chromosomes to events one chromosome at a time
        emap = zeros(Nevents,Nsamples); % allocate storage for expanded event map
        for c = 1:Nchr
            emap(chrns_e{c},:) = double(E.dat(chrns_e{c},idx_mat(c,:)));
        end
        % add up co-occurences with matrix multiplication
        cooc = emap * emap';

        %%%%% begin method fork
        % compare with observed, count equality and each direction of inequality
        le_counts = le_counts + (cooc <= obs_tot);
        ge_counts = ge_counts + (cooc >= obs_tot);
        eq_counts = eq_counts + (cooc == obs_tot);
        % accumulate le & ge counts for lineage subsets
        for l = 1:Nlineages
            submap = emap(:,lindices{l});
            lcooc = submap * submap';
            le_counts_byclass(:,:,l) = le_counts_byclass(:,:,l) + (lcooc <= lin_obs(:,:,l));
            ge_counts_byclass(:,:,l) = ge_counts_byclass(:,:,l) + (lcooc >= lin_obs(:,:,l));
            eq_counts_byclass(:,:,l) = eq_counts_byclass(:,:,l) + (lcooc == lin_obs(:,:,l));
        end
        %%%%% end method fork

    end
    toc
end % loop over chunk files

%% calculate statistics for pairs of events

% calculate number of pairs on different chromosomes
Npairs = 0;
for c = 1:Nchr
    Npairs = Npairs + sum(E.event.chrn(c) < E.event.chrn);
end

pscnt = options.pcount;

% allocate storage for lineage-specific p-values
pairs_idx = zeros(Npairs,2);            % list of unique event pair indices

%% count events for lineage-specfic correlation and anti-correlation p-values
if Nlineages > 0
    p_list_corr = zeros(Npairs,Nlineages); % correlation p-values per lineage
    p_list_anti = zeros(Npairs,Nlineages); % anti-correlation p-values per lineage

    %%%%% begin method fork
    % eg/gt/lt has nothing to do here
    %%%%% end method fork

    verbose('calculating lineage-specific p-values',10);
    tic
    for l = 1:Nlineages
        s = 1; % pair index
        p_list_cpow = zeros(Npairs,1);  % min fisher exact p-value for correlation
        p_list_apow = zeros(Npairs,1);  % min fisher exact p-value for anti-correlation
        for i = 1:Nevents
            for j = i+1:Nevents
                if E.event.chrn(i) ~= E.event.chrn(j)
                    % maximum power calculation
                    lx = lindices{l}; % indices of samples in this lineage
                    [p_list_cpow(s),p_list_apow(s)] = max_fish_power(length(lx), sum(E.dat(i,lx)), sum(E.dat(j,lx))); 
                    % calculate p-values for both correlation and anti-correlation
                    if options.split_eq
                        % historical
                        p_list_corr(s,l) = (ge_counts_byclass(i,j,l) - eq_counts_byclass(i,j,l)/2 + pscnt) / (Nperms + pscnt); 
                        p_list_anti(s,l) = (le_counts_byclass(i,j,l) - eq_counts_byclass(i,j,l)/2 + pscnt) / (Nperms + pscnt);
                    else
                        % the correct way of dealing with permuted==observed
                        p_list_corr(s,l) = (ge_counts_byclass(i,j,l) + pscnt) / (Nperms + pscnt);
                        p_list_anti(s,l) = (le_counts_byclass(i,j,l) + pscnt) / (Nperms + pscnt);
                    end
                    pairs_idx(s,:) = [i,j];
                    s = s + 1;
                end
            end
        end

        % output selected lineage tables
        output_pair_p(E.event,[p_list_anti(:,l),pairs_idx,p_list_apow],...
                      fullfile(save_dir,['anticorr_pair.',linames{l},ext,'.txt']),...
                      1,options.sig_thresh,options.power_thresh);
        output_pair_p(E.event,[p_list_corr(:,l),pairs_idx,p_list_cpow],...
                      fullfile(save_dir,['correlate_pair.',linames{l},ext,'.txt']),...
                      1,options.sig_thresh,options.power_thresh);
    end
    toc
end

%% count events for overall correlation and anti-correlation p-values

% sum permuted per-lineage event coincidences into overall permuted event coincidences
verbose('counting overall event co-occurrences',10)

% allocate storage for lineage-specific p-values
pairs_idx = zeros(Npairs,2);   % peak indices
p_corr = zeros(Npairs,1);
p_anti = zeros(Npairs,1);
p_cpow = zeros(Npairs,1);  % min fisher exact p-value for correlation
p_apow = zeros(Npairs,1);  % min fisher exact p-value for anti-correlation

% calculate p-values

%%%%% begin method fork
% eq/lt/gt method has nothing to do
%%%%% end method fork

verbose('calculating overall p-values for co-occurrences',10);
tic
s = 1; % pair index
for i = 1:Nevents-1
    for j = i+1:Nevents
        if E.event.chrn(i) ~= E.event.chrn(j)
            % maximum power calculation
            [p_cpow(s),p_apow(s)] = max_fish_power(Nsamples,sum(E.dat(i,:)),sum(E.dat(j,:)));
            % p value calculations
            if options.split_eq
                p_corr(s) = (ge_counts(i,j) - eq_counts(i,j)/2 + pscnt) / (Nperms + pscnt);
                p_anti(s) = (le_counts(i,j) - eq_counts(i,j)/2 + pscnt) / (Nperms + pscnt); 
            else
                p_corr(s) = (ge_counts(i,j) + pscnt) / (Nperms + pscnt);
                p_anti(s) = (le_counts(i,j) + pscnt) / (Nperms + pscnt);
            end
            % save pair index 
            pairs_idx(s,:) = [i,j];
            s = s+1;
        end
    end
end
toc

%% save results files
verbose('saving results',10);toc

% save binaries for forensics
save(fullfile(save_dir,['pair_results',ext,'.mat']),'pairs_idx','p_corr','p_anti','p_cpow','p_apow');
%!save(fullfile(save_dir,['pair_perm_tot',ext,'.mat'],'perm_tot');
save(fullfile(save_dir,['pair_obs_tot',ext,'.mat']),'obs_tot');

% save overall significance results files
output_pair_p(E.event,[p_anti,pairs_idx,p_apow],...
                  fullfile(save_dir,['anticorr_pair',ext,'.txt']),...
                  1,options.sig_thresh,options.power_thresh);
output_pair_p(E.event,[p_corr,pairs_idx,p_cpow],...
                  fullfile(save_dir,['correlate_pair',ext,'.txt']),...
                  1,options.sig_thresh,options.power_thresh);

% save function WS for debugging and downstream analyses
save(fullfile(save_dir,['analyze_pairs2_ws',ext,'.mat']));


end % function

