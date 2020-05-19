function [D,margs_sort,new_samples] = corrperm_prep(D,regs,work_dir,options)
% prepare inputs for permutation tests: module level function
%
%   create lineage "bins"
%   call events
%   save outputs to work_dir
%
% derived from merger of functions 
% /xchip/gistic/Travis/Correlations/LSF2/0816/chrm_metro_restart_lsf2_lineage
% /xchip/gistic/Travis/Correlations/LSF2/0816/chrm_metro_restart_lsf2_lineage_hl

%% process optional parameters
if ~exist('options','var') || isempty(options)
    options = struct;
end
options = impose_default_value(options,'permclass_sisfield',[]); % no permutation classes

options = impose_default_value(options,'hilevel',false);
options = impose_default_value(options,'minclass_samples',17); % 15 for high-level
options = impose_default_value(options,'t_amp',0.3);
options = impose_default_value(options,'t_del',0.3);
options = impose_default_value(options,'broad_len_cutoff',0.50);
options = impose_default_value(options,'event_thresh',0);
options = impose_default_value(options,'max_disruption',[Inf,Inf]);

% get number of chromosomes from D
Nchr = max(D.chrn);

% create working directory if necessary
if ~exist(work_dir,'dir')
    mkdir(work_dir)
end

%% calculate disruption score for each chromosome/sample/SCNA type
verbose('scoring observed marginals',20)

% set parameters for reconstruct_genomes
params = struct;
params.use_segarray = true;
params.broad_or_focal = 'focal';
params.broad_len_cutoff = options.broad_len_cutoff;
params.t_amp = options.t_amp;
params.t_del = options.t_del;
params.column_to_add = 4; %!!! move into function

if options.hilevel
    % high-level: focal amps, all-event dels
    focal_genome = reconstruct_genomes(D.Qs,params);
    F.dat = focal_genome.amp-focal_genome.del;
    scores1.amp = F.dat>(params.t_amp+options.event_thresh);
    params.broad_len_cutoff = 2.1;
    broad_genome = reconstruct_genomes(D.Qs,params);
    B.dat = broad_genome.amp-broad_genome.del;
    scores1.del = B.dat<-(params.t_del+options.event_thresh);
else
    % low level: focal amps and dels
    scores1 = reconstruct_genomes(D.Qs,params);
end

% calculate marginal disruption score
Nsamples = size(D.dat,2);
margs = zeros(Nchr,Nsamples,2);
Dchrn = SegArray(D.chrn);
fields = {'amp','del'};
for i = 1:Nchr
    idx = Dchrn==i;
    for k = 1:2
        % score by summing events
        margs(i,:,k) = sum(scores1.(fields{k})(idx,:)~=0);%!!! ./sum(idx);
        %!!! TODO have normalized-by-sum code variation to allow equal
        %!!!      contributions for sample regardless of disruption
    end
end

%% eliminate samples with too much disruption
good_idx = find(sum(margs(:,:,2))<options.max_disruption(2) & ...
                sum(margs(:,:,1))<options.max_disruption(1) );
D = reorder_D_cols(D,good_idx);
margs = margs(:,good_idx,:);
Nsamples = size(D.dat,2);

%% form permutation classes (lineages) and eliminate unpowered lineages
[D,samples,reorder] = create_D_sample_bins(D,options.minclass_samples,options.permclass_sisfield);
margs = margs(:,reorder,:);
marg = squeeze(sum(margs));
% stash amp/del marginals in D
D.amp_disrupt = marg(:,1)';
D.del_disrupt = marg(:,2)';
add_D_field(D,'sample',{'amp_disrupt','del_disrupt'});

%% put the observed samples in optimal order for fitting during permutations
order = 1:Nsamples;
new_samples = samples;  % init reordered permutation class indices
% loop over classes
for i = 1:length(samples)
    [~,I] = sort(sum(marg(samples{i},:),2));
    new_samples{i} = sort(samples{i}(I)); %!
%!  new_samples{i} = samples{i}(I);
    order(samples{i}) = order(samples{i}(I));
end
% sort the marginals
margs_sort = margs(:,order,:);
% sort the copy number data
D = reorder_D_cols(D,order);

%% call events at peaks across samples
% this is not strictly necessary here - event scoring can be done
% at the analysis stage after permutations

% get peak center snp indices
amps = [regs{1}.peak];
dels = [regs{2}.peak]; 

params.broad_len_cutoff = options.broad_len_cutoff;
[ Binary,Binary_amps,Binary_dels ] = score_cooccurance(D.Qs,amps,dels,params,...
                                            options.hilevel,options.event_thresh);

%% save input files to the permutation/evaluation stages in working directory
verbose('saving permutation input files.',20)

% inputs to permutation engine
save(fullfile(work_dir,'new_samples.mat'),'new_samples');
save(fullfile(work_dir,'margs.mat'),'margs_sort');

%! permutation input as single package
H = struct;
H.sdesc = D.sdesc;
H.margs = margs_sort;
H.classes = new_samples;
save(fullfile(work_dir,'H.mat'),'H');

% inputs to evaluator
save(fullfile(work_dir,'Binary_amps.mat'),'Binary_amps');
save(fullfile(work_dir,'Binary_dels.mat'),'Binary_dels');
save(fullfile(work_dir,'D.mat'),'D');
save(fullfile(work_dir,'peak_regs.mat'),'regs');

%! event map as a single package
E = struct;
E.sdesc = D.sdesc;
E.event = struct;
E.event.name = [strcat('amp_',{regs{1}.name}),strcat('del_',{regs{1}.name})]';
E.event.type = repmat(1,length(regs{1}),1);repmat(2,length(regs{2}),1);
E.event.chrn = [regs{1}.chrn,regs{2}.chrn]';
E.dat = [Binary_amps;Binary_dels];
save(fullfile(work_dir,'Emap.mat'),'E');
