function [idx_cell,stat_finals,stats] = ampdel_tempering(H,iters,opts,varargin)
%CORRPERM_AMPDEL_TEMPERING iterate over rounds of simulated tempering
%
%   [IDX_CELL,STAT_FINALS,STATS] = ampdel_tempering(H,ITERS,OPTS,...)
%
% H is the disruption profile structure of the data. ITERS is the number of trials.
% OPTS are permutation options.
%
% returns three cell arrays with an element for each permutation:
% IDX_CELL contains permutation indices for each run (integers, Nchr X Nsamples)
% STAT_FINALS contains containing error vectors {amp|del} for each trial
% STATS contains structs of loosely defined statistics for analyzing the runs


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
[Nchr,Nsamples,Nalt] = size(H.margs);

% find minimum non-zero disruption values for each chromosome for amps and dels
min_bin = nan(Nchr,2);
for i = 1:Nchr      % loop over chromosomes
    for j = 1:Nalt  % amp/del
        test = H.margs(i,:,j);
        test(test==0) = Inf;
        min_bin(i,j) = min(test) + opts.minbinmin;
    end
end

rand_margs_cell = cell(1,iters);
stats = cell(1,iters);
stat_finals = cell(1,iters);
idx_cell = cell(1,iters);

margs1 = zeros(size(H.margs));
idx_mat0 = repmat(1:Nsamples,Nchr,1);
idx_mat1 = repmat(1:Nsamples,Nchr,1);

%% iterate over specified number of permutations
tic
for d = 1:iters
    % randomly permute sample labels within permutation classes
    for k = 1:length(H.pcx)   % loop over permutation classses
        for j = 1:Nchr              % loop over chromosomes
            r = randperm(length(H.pcx{k}));
            margs1(j,H.pcx{k},:) = H.margs(j,H.pcx{k}(r),:);
            idx_mat1(j,H.pcx{k}) = idx_mat0(j,H.pcx{k}(r));
        end
    end
    
    % sort permuted samples by disruption within each class to optimize fitting
    rand_marg = squeeze(sum(sum(margs1,1),3)); % randomized sample marginals
    rand_margs = NaN(size(margs1));
    idx_mat_s = zeros(size(idx_mat1));
    for i = 1:length(H.pcx)
        [~,I1] = sort(rand_marg(H.pcx{i}));
        rand_margs(:,H.pcx{i},:) = margs1(:,H.pcx{i}(I1),:);
        idx_mat_s(:,H.pcx{i}) = idx_mat1(:,H.pcx{i}(I1));
    end
    
    % run the annealing schedule
    [idx_mat2,rand_margs,stat_final,stat] = corrperm_tempering_schedule(...
                        H.margs,rand_margs,idx_mat_s,[10,10],H.pcx,min_bin,opts);

    rand_margs_cell{d} = rand_margs ;
    stats{d} = stat;
    stat_finals{d} = stat_final;
    idx_cell{d} = idx_mat2;
    verbose('iteration %d: %0.1f seconds',20,d,toc);

end % permutation iteration
total_time = toc;
verbose('%d permutations in %0.1f seconds: %0.1f permutations/hour\n',10,...
                        iters,total_time,iters*3600/total_time);
verbose('- exiting permutation run -',10);
end
