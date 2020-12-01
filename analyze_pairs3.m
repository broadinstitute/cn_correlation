function analyze_pairs3(E,perm_dir,save_dir,options)
% memory-conserving version of corrperm_analyze_pairs that accumulates a full histogram of
% permuted event counts for each pair of events
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
options = impose_default_value(options,'perm_file_mask','idx_cell*');% set to split obs==permuted between tail and body 
options = impose_default_value(options,'lineage_out',{});
options = impose_default_value(options,'analyze_lineages',false);
options = impose_default_value(options,'lineage_out',E.pcname);    % list of lineage subsets to analyze

options

% input warnings
%!if options.split_eq
%!     warning('The ''split_eq'' way of calculating P-values is deprecated.');
%!end
if options.pcount == 0
     warning('Not using a pseudocount is deprecated.');
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

% initialize lineage breakdown working variables
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
% This version of the analysis aggregates a histogram of co-oocurrence counts
% for each each event pair.
%

% overall
% use Nevents x Nevents array of histograms
maxcount = max(sum(E.dat,2));  % maximum possible co-occurrence
cochist = zeros(Nevents,Nevents,maxcount+1);
% by lineage
if Nlineages > 0
    cochist_byclass = cell(Nlineages,1);
    for l = 1:Nlineages
        maxcount = max(sum(E.dat(:,lindices{l}),2));
        cochist_byclass{l} = zeros(Nevents,Nevents,maxcount+1);
    end
end
%%%%%
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
        %% map chromosome-swap indices to event swaps
        % allocate storage
        emap = zeros(Nevents,Nsamples); % allocate storage for expanded event map
        % map events one chromosome at a time
        for c = 1:Nchr
            emap(chrns_e{c},:) = double(E.dat(chrns_e{c},idx_mat(c,:)));
        end
        % add up co-occurences with matrix multiplication
        cooc = emap * emap';
        %%%%% begin method fork
        % update histograms
        cochist = count_binsum(cochist,cooc);

        % add up lineage subsets of co-occurences and update histograms
        for l = 1:Nlineages
            submap = emap(:,lindices{l});
            lcooc = submap * submap';
            cochist_byclass{l} = count_binsum(cochist_byclass{l},lcooc);
        end
        %%%%% end method fork
    end
    toc
end % loop over chunk files

%% calculate statistics for pairs of events

% calculate number of pairs on different chromosomes
Npairs = 0;
for i = 1:Nevents
    Npairs = Npairs + sum(E.event.chrn(i) < E.event.chrn);
end

pscnt = options.pcount;

% allocate storage for lineage-specific p-values
pairs_idx = zeros(Npairs,2);            % list of unique event pair indices

%% count events for lineage-specfic correlation and anti-correlation p-values
if Nlineages > 0
    p_list_corr = zeros(Npairs,Nlineages); % correlation p-values per lineage
    p_list_anti = zeros(Npairs,Nlineages); % anti-correlation p-values per lineage
    
    %%%%% begin method fork

    % calculate by-lineage ge/le/eq values
    le_counts_byclass = zeros(Nevents,Nevents,Nlineages);
    ge_counts_byclass = zeros(Nevents,Nevents,Nlineages);
    eq_counts_byclass = zeros(Nevents,Nevents,Nlineages);

    verbose('calculating lineage-specific le/ge/eq counts',10);
    tic
    for l = 1:Nlineages
       for i = 1:Nevents
            for j = i+1:Nevents
                this_obs = lin_obs(i,j,l)+1;
                this_perm = cochist_byclass{l}(i,j,:);
                eq_counts_byclass(i,j,l) = this_perm(this_obs);
                le_counts_byclass(i,j,l) = sum(this_perm(1:this_obs));
                ge_counts_byclass(i,j,l) = sum(this_perm(this_obs:end));
            end
       end
    end
    toc
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

        % output correlation/anticorrelation tables for selected lineages
        output_pair_p(E.event,[p_list_anti(:,l),pairs_idx,p_list_apow],...
                      fullfile(save_dir,['anticorr_pair.',linames{l},ext,'.txt']),...
                      1,options.sig_thresh,options.power_thresh);
        output_pair_p(E.event,[p_list_corr(:,l),pairs_idx,p_list_cpow],...
                      fullfile(save_dir,['correlate_pair.',linames{l},ext,'.txt']),...
                      1,options.sig_thresh,options.power_thresh);
    end % 1:Nlineages
    toc
end % if Nlineages > 0

%% count events for overall correlation and anti-correlation p-values

% allocate storage for lineage-specific p-values
pairs_idx = zeros(Npairs,2);   % peak indices
p_corr = zeros(Npairs,1);
p_anti = zeros(Npairs,1);
p_cpow = zeros(Npairs,1);  % min fisher exact p-value for correlation
p_apow = zeros(Npairs,1);  % min fisher exact p-value for anti-correlation
% calculate p-values

%%%%% begin method fork

%% calculate counts for p-values
verbose('calculating overall le/ge/eq counts',10);
tic

le_counts = zeros(Nevents,Nevents);
ge_counts = zeros(Nevents,Nevents);
eq_counts = zeros(Nevents,Nevents);

for i = 1:Nevents-1
    for j = i+1:Nevents
        this_obs = obs_tot(i,j)+1;
        this_perm = cochist(i,j,:);
        eq_counts(i,j) = this_perm(this_obs);
        le_counts(i,j) = sum(this_perm(1:this_obs));
        ge_counts(i,j) = sum(this_perm(this_obs:end));
    end
end
toc

%%%%% end method fork

verbose('calculating overall co-occurrences p-values',10);
tic
s = 1; % pair index
for i = 1:Nevents-1
    for j = i+1:Nevents
        if E.event.chrn(i) ~= E.event.chrn(j)
            %% maximum power calculation
            [p_cpow(s),p_apow(s)] = max_fish_power(Nsamples,sum(E.dat(i,:)),sum(E.dat(j,:)));
            %% p value calculations
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
verbose('saving results',10)

% save binaries for forensics
save(fullfile(save_dir,['pair_results',ext,'.mat']),'pairs_idx','p_corr','p_anti','p_cpow','p_apow');
save(fullfile(save_dir,['cochist',ext,'.mat']),'cochist');
save(fullfile(save_dir,['pair_obs_tot',ext,'.mat']),'obs_tot');

% save overall significance results files
output_pair_p(E.event,[p_anti,pairs_idx,p_apow],...
                  fullfile(save_dir,['anticorr_pair',ext,'.txt']),...
                  1,options.sig_thresh,options.power_thresh);
output_pair_p(E.event,[p_corr,pairs_idx,p_cpow],...
                  fullfile(save_dir,['correlate_pair',ext,'.txt']),...
                  1,options.sig_thresh,options.power_thresh);

% save function WS for debugging and downstream analyses
save(fullfile(save_dir,['analyze_pairs3_ws',ext,'.mat']));

end % main function


%% subfunction for updating histogram

function cochist = count_binsum(cochist,cocounts)
% cochist is histogram of counts across permutations (Nevents x Nevents x max_counts)
% cocounts are the coocurrences for a single iteration (Nevents x Nevents)

    % if count excededs size of histogram, reallocate (shouldn't happen, but does)
    if max(max(cocounts)) >= size(cochist,3)
        cochist(1,1,max(max(cocounts))+1) = 0;
    end

    % bump histogram based on co-occurences
    Nevents = size(cocounts,1);
    for i = 1:Nevents-1
        for j = i+1:Nevents
            c = cocounts(i,j);
            cochist(i,j,c+1) = cochist(i,j,c+1)+1;
        end
    end
end % subfunction
