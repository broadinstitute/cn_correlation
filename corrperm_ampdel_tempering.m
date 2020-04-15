function [rand_margs_cell,idx_cell,stat_finals,stats] = corrperm_ampdel_tempering(margs,samples,iters,opts)
%CORRPERM_AMPDEL_TEMPERING iterate over rounds of simulated tempered annealing
%
%   [RAND_MARGS_CELL,IDX_CELL,STATS] = 
%            corrperm_ampdel_tempering(MARGS,IDX_MAT,SAMPLES,ITERS,OPTS)
%
% MARGS are the observed marginals for each chromosome of each sample for
% every SCNA type (Nchr X Nsamples X amp|del). IDX_MAT is a set of input
% permutation indices.

% returns three cell arrays with an element for each permutation:
% RAND_MARGS_CELL are the permutation marginals (each element disruption sum, 
% Nchr X Nsamples X amp|del); IDX_CELL are the permutaion indices (integers, 
% Nchr X Nsamples); STATS are statistics (loosely defined).


% derived from:
% /xchip/gistic/Travis/Correlations/LSF2/0816/annealing_permutations_lsf2
% /xchip/gistic/Travis/Correlations/LSF2/0816/annealing_permutations_lsf2_hl

% optional parameter processing
if ~exist('opts','var') || isempty(opts)
    opts = struct;
end
opts = impose_default_value(opts,'binminbin',0);

verbose('- initiating permutation run of %d trials -',10,iters);

% number of chromosomes
[Nchr,Nsamples,Nalt] = size(margs);

% find minimum non-zero disruption values for each chromosome for amps and dels
min_bin = nan(Nchr,2);
for i = 1:Nchr      % loop over chromosomes
    for j = 1:Nalt  % amp/del
        test = margs(i,:,j);
        test(test==0) = Inf;
        min_bin(i,j) = min(test) + opts.minbinmin;
    end
end

% randomize random number generator for this chunk
if str2num(regexprep(version,'\.[0-9]+\.[0-9]+ \(R.+\)$','')) >= 8.0
    rng('shuffle') % Matlab R2012b
else
    rand('seed',rem(now*10e9,1e6));disp(randperm(10)); %! Matlab R2010b: no 'rng', seed rng from wall clock
%!    rand('seed',now); %!!! does not effectively randomize the seed
end
%min_bin = 100.*ones(23,2);%min_bin;

rand_margs_cell = cell(1,iters);
stats = cell(1,iters);
stat_finals = cell(1,iters);
idx_cell = cell(1,iters);

margs1 = zeros(size(margs));
idx_mat0 = repmat(1:Nsamples,Nchr,1);
idx_mat1 = repmat(1:Nsamples,Nchr,1);

%% iterate over specified number of permutations
%matlabpool open 4
%parfor d = 1:iters
tic
for d = 1:iters
    % randomly permute sample labels within permutation classes
    for k = 1:length(samples)   % loop over permutation classses
        for j = 1:Nchr              % loop over chromosomes  !!!HUMGEN
            r = randperm(length(samples{k}));
            margs1(j,samples{k},:) = margs(j,samples{k}(r),:);
            idx_mat1(j,samples{k}) = idx_mat0(j,samples{k}(r));
        end
    end
    
    % sort permuted samples by disruption within each class to optimize fitting
    rand_marg = squeeze(sum(sum(margs1,1),3)); % randomized sample marginals
    rand_margs = NaN(size(margs1));
    idx_mat_s = zeros(size(idx_mat1));
    for i = 1:length(samples)
        [~,I1] = sort(rand_marg(samples{i}));
        rand_margs(:,samples{i},:) = margs1(:,samples{i}(I1),:);
        idx_mat_s(:,samples{i}) = idx_mat1(:,samples{i}(I1));
    end
    
    % run the annealing schedule
    [idx_mat2,rand_margs,stat_final,stat] = corrperm_tempering_schedule(...
                        margs,rand_margs,idx_mat_s,[10,10],samples,min_bin,opts);

    rand_margs_cell{d} = rand_margs ;
    stats{d} = stat;
    stat_finals{d} = stat_final;
    idx_cell{d} = idx_mat2;
    verbose('iteration %d: %0.1f seconds',20,d,toc);

end % permutation iteration
total_time = toc;
verbose('%d permutations in %0.1f seconds: %0.1f permutations/hour\n',10,...
                        iters,total_time,iters*3600/total_time);
%matlabpool close
verbose('- exiting permutation run -',10);


