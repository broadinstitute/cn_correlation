function [idx_mat2,rand_margs2,stat_final,stat_data] = corrperm_tempering_schedule(margs,rand_margs,idx_mat1,step,samples,min_bin,opts)
%CORRPERM_TEMPERING_SCHEDULE - parameterized basic simulated tempering
%
% [IDX_MAT2,RAND_MARGS2,STAT_FINAL,STAT_DATA] = corrperm_tempering_schedule(...
%           MARGS,RAND_MARGS,IDX_MAT1,STEP,SAMPLES,MIN_BIN,OPTS)
%
%
% Use simulated tempering to exchange chromosomes of samples while
% controlling for disruption. The indices in IDX_MAT1 and marginals in
% RAND_MARGS are transformed into IDX_MAT2 and RAND_MARGS2.
%
% TO DO finish documentation

% derived from:
%   /xchip/gistic/Travis/Correlations/LSF2/cooldown_cycle_ratio_lsf2.m
%   /xchip/gistic/Travis/Correlations/LSF2/cooldown_cycle_ratio_lsf2_hl.m

if ~strcmp(opts.algorithm,'basic_tempered_annealing')
    error('algorithm mismatch');
end

stat_data = [];

%% initial cooldown
[rand_margs2,idx_mat2,stat_final,new_stats] = corrperm_cooldown_ratio(margs,rand_margs,...
                        idx_mat1,opts.initial.iters,opts.initial.temp,...
                        opts.initial.step,samples,min_bin,opts);
stat_data = catstats(stat_data,new_stats);

%% tempering cycle(s)
if isfield(opts,'cycle')
    for i = 1:length(opts.cycle)
        test1 = stat_final(end,:);
        % adjust start temperature to "melt" the higher of amps or dels
        if test1(1)>test1(2)
            starts = opts.cycle(i).heatup{1};
        else
            starts = opts.cycle(i).heatup{2};
        end

        [rand_margs2,idx_mat2,stat_final,new_stats] = corrperm_cooldown_ratio(...
                            margs,rand_margs2,idx_mat2,opts.cycle(i).iters,...
                            starts,opts.cycle(i).step,samples,min_bin,opts);
        stat_data = catstats(stat_data,new_stats);
    end
end

%% final cooldown
[rand_margs2,idx_mat2,stat_final,new_stats] = corrperm_cooldown_ratio(...
                        margs,rand_margs2,idx_mat2,opts.final.iters,...
                        opts.final.temp,opts.final.step,samples,min_bin,opts);
stat_data = catstats(stat_data,new_stats);

end % function

%% subfunction: catenate statistics
% updates a struct of vectors
function data = catstats(data,new_data)
    if isempty(data)
        data = new_data;
    else
        for fld = fields(new_data)'
            f = fld{1};
            data.(f) = [data.(f);new_data.(f)];
        end
    end
end

