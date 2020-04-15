function [idx_mat2,rand_margs2,stat_final,stat_data] = corrperm_lockstep_schedule(margs,rand_margs,idx_mat1,step,samples,min_bin,opts)
%CORRPERM_TEMPERING_SCHEDULE - parameterized basic simulated tempering
%
% [IDX_MAT2,RAND_MARGS2,STAT_FINAL,STAT_DATA] = corrperm_tempering_schedule(...
%           MARGS,RAND_MARGS,IDX_MAT1,STEP,SAMPLES,MIN_BIN,OPTS)
%
%

% derived from:
%   /xchip/gistic/Travis/Correlations/LSF2/cooldown_cycle_ratio_lsf2.m
%   /xchip/gistic/Travis/Correlations/LSF2/cooldown_cycle_ratio_lsf2_hl.m

if ~strcmp(opts.algorithm,'basic_lockstep_annealing')
    error('algorithm mismatch');
end

stat_data = [];

%% initial cooldown
if isfield(opts,'initial') && opts.initial.iters > 0
    [rand_margs2,idx_mat2,stat_final,new_stats] = corrperm_cooldown_ratio(...
                            margs,rand_margs,idx_mat1,opts.initial.iters,...
                            opts.initial.temp,opts.initial.step,samples,min_bin,opts);
    stat_data = catstats(stat_data,new_stats);
end
temps = stat_data.temp(end,:);
%% tempering cycle(s)
if isfield(opts,'cycle')
    for i = 1:length(opts.cycle)
        tadjust = opts.cycle(i).coolfactor;
        if isscalar(tadjust)
            tadjust = [tadjust,1];
        end
        test1 = stat_final(end,:);
        % adjust start temperature to "melt" the higher of amps or dels
        if test1(1)>test1(2)
            if temps(1) < temps(2)
                temps = [1 1] .* prod(temps)^(1/2);
            end  
            temps = temps .* tadjust;               
        else
            if temps(2) < temps(1)
                temps = [1 1] .* prod(temps)^(1/2);
            end  
            temps = temps .* fliplr(tadjust);
        end

        [rand_margs2,idx_mat2,stat_final,new_stats] = corrperm_cooldown_ratio(...
                            margs,rand_margs2,idx_mat2,opts.cycle(i).iters,...
                            temps,opts.cycle(i).step,samples,min_bin,opts);
        stat_data = catstats(stat_data,new_stats);
    end
end

%% final cooldown
if isfield(opts,'final') && opts.final.iters > 0
    [rand_margs2,idx_mat2,stat_final,new_stats] = corrperm_cooldown_ratio(...
                            margs,rand_margs2,idx_mat2,opts.final.iters,...
                            opts.final.temp,opts.final.step,samples,min_bin,opts);
    stat_data = catstats(stat_data,new_stats);
end

%% subfunction: catenate statistics
function data = catstats(data,new_data)

if isempty(data)
    data = new_data;
else
    for fld = fields(new_data)'
        f = fld{1};
        data.(f) = [data.(f);new_data.(f)];
    end
end
