function [rand_margs,idx_mat,stat_final,data] = corrperm_cooldown_ratio(...
                margs,rand_margs,idx_mat,g,Temp,step,samples,min_bin,opts)
%COOLDOWN_RATIO do a single cooling phase for simulated annealing/tempering
%
%[RAND_MARGS,IDX_MAT,STAT_FINAL,DATA] = corrperm_cooldown_ratio( ...
%                                       MARGS,RAND_MARGS,IDX_MAT,G, ...
%                                       TEMP,STEP,SAMPLES,MIN_BIN,OPTS)
%
% MARGS - observed marginals by chromosome, chrom X sample X amp|del
% RAND_MARGS - permuted marginals by chromosome, chrom X sample X amp|del
% IDX_MAT - scrambled permutation index, chrom X sample
% G - number of iterations of a permutation class to perform at each temperature
%     (effectively, the length of time to cool for)
% TEMP - start temperature NOTE: cooler temperature is higher. more like energy ~ exp(1/RT) 
% STEP - factor to mutiply temperature by when adjusting it
% SAMPLES - definition of permutation classes (cell array of index vectors)
% MIN_BIN - minimum non-zero disruption value, chromosome X amp|del
% OPTS - options
%           OPTS.stat_interval how often to gather data about the run
%           (default 1000)
%           OPTS.class_subiters how many swaps to do in the same class per
%           iteration defined by G
% return:
%   RAND_MARGS new permuted marginals
%   IDX_MAT new permuted index
%   STAT_FINAL is an amp/del vector of the final absolute value of the error
%       between observed and annealed data
%   DATA is a struct containing useful information about the run sampled at
%       times determined by the OPTS.stat_interval parameter:
%       



% derived from:
% /xchip/gistic/Travis/Correlations/LSF2/cooldown_ratio_lsf2.m

if ~exist('opts','var')
    opts = struct;
end

opts = impose_default_value(opts,'stat_interval',1000); % interval for recording statistics
opts = impose_default_value(opts,'class_subiters',100); % number of iterations in same class

Temp_resets = 1; %!!! updated but unused, DELETE ME

% sum chromosomes into marg (samples X amp/del)
marg = squeeze(sum(margs));             % observed
rand_marg = squeeze(sum(rand_margs));   % permuted

% initialize data structure for statstics
Nstats = floor(g/opts.stat_interval);
data = struct('temp',nan(Nstats,2),...
              'stats',nan(Nstats,2),...
              'msd',nan(Nstats,2),...
              'ratios',nan(Nstats,1),...
              'backsteps',nan(Nstats,1));
% initialize statistics
pp = 1; % index for stats
%{
data.stats(pp,:) = sum(((rand_marg-marg)./max(marg,2500)).^2);
data.temp(pp,:) = Temp;
data.ratios(pp) = 0;
data.backsteps(pp) = 0;
%}
Nchr = size(margs,1);


%}

%!UNUSED stat_initial = sum((rand_marg_initial-marg).^2);
%!UNNCESSARY rand_marg = rand_marg_initial;
%!UNUSED u = 0;
k = 0;
a = 0;
cellsz = cellfun(@length,samples,'uni',false);
cellsz = cell2mat(cellsz);

%% iterate over temperature steps 
for l = 1:g
    % pre-select (for efficiency) opts.class_subiters sample pairs from the same chromosome
    % and the same permutation class (lineage) as proposed swaps
    samps = randsample(length(samples),1,true,cellsz); % randomlypick a class
    n = length(samples{samps}); % number of samples in class
    ss = samples{samps}(randi(n,1,opts.class_subiters));
    tt = samples{samps}(randi(n,1,opts.class_subiters));
    cc = randi(Nchr,1,opts.class_subiters);
    
    % opts.class_subiters at same (temperature,permutation class)
    for j = 1:opts.class_subiters
        % index preselected pair/chromosome
        s = ss(j);
        t = tt(j);
        c = cc(j);
        
        %% evaluate effect of swap for amplifications (w/o actually swapping)
        E_amp = 0;
        delta_amp = rand_margs(c,s,1)-rand_margs(c,t,1);
        if delta_amp ~=0
            % denominators for normalizing disruption
            denorm_as = max(marg(s,1),min_bin(c,1));
            denorm_at = max(marg(t,1),min_bin(c,1));
            if marg(s,1) ~= marg(t,1) || denorm_as ~= denorm_at
                % marginal errors
                margerr_as = rand_marg(s,1)-marg(s,1);
                margerr_at = rand_marg(t,1)-marg(t,1);
                % disruption values, old & new
                old_at_value = margerr_at ./ denorm_at;
                new_at_value = (margerr_at+delta_amp) ./ denorm_at;
                old_as_value = margerr_as ./ denorm_as;
                new_as_value = (margerr_as-delta_amp) ./ denorm_as;
                % compute amplification contribution to the energy
                E_amp = -abs(old_at_value)+abs(new_at_value)-abs(old_as_value)+abs(new_as_value);
            end
        end
        
        %% evaluate effect of swap for deletions
        E_del = 0;
        delta_del = rand_margs(c,s,2)-rand_margs(c,t,2);
        if delta_del ~=0
            % denominators for normalizing disruption
            denorm_ds = max(marg(s,2),min_bin(c,1));
            denorm_dt = max(marg(t,2),min_bin(c,1));
            if marg(s,2) ~= marg(t,2) || denorm_ds ~= denorm_dt
                % marginal errors
                margerr_ds = rand_marg(s,2)-marg(s,2);
                margerr_dt = rand_marg(t,2)-marg(t,2);
                % disruption values, old & new
                old_dt_value = margerr_dt ./ denorm_dt;
                new_dt_value = (margerr_dt+delta_del) ./ denorm_dt;
                old_ds_value = margerr_ds ./ denorm_ds;
                new_ds_value = (margerr_ds-delta_del)./denorm_ds;
                % compute amp energy
                E_del = -abs(old_dt_value)+abs(new_dt_value)-abs(old_ds_value)+abs(new_ds_value);
            end
        end
        
        %% accept or reject the proposed swap based on total energy
        E = Temp(1)*E_amp + Temp(2)*E_del;

        rn = rand(1);
        if E <= 1e-5
            %!!! temporary test/breakpoint location
%!          if E > 0
%!              fprintf('SNORT!!! E =%d\n',E);
%!          end
            rand_marg(s,1) = rand_marg(s,1)-delta_amp;
            rand_marg(t,1) = rand_marg(t,1)+delta_amp;
            rand_marg(s,2) = rand_marg(s,2)-delta_del;
            rand_marg(t,2) = rand_marg(t,2)+delta_del;

            rand_margs(c,s,1) = rand_margs(c,s,1)-delta_amp;
            rand_margs(c,s,2) = rand_margs(c,s,2)-delta_del;
            rand_margs(c,t,1) = rand_margs(c,t,1)+delta_amp;
            rand_margs(c,t,2) = rand_margs(c,t,2)+delta_del;
            temp = idx_mat(c,s);
            idx_mat(c,s) = idx_mat(c,t);
            idx_mat(c,t) = temp;
            k = k+1;
        elseif E < rn
            rand_marg(s,1) = rand_marg(s,1)-delta_amp;
            rand_marg(t,1) = rand_marg(t,1)+delta_amp;
            rand_marg(s,2) = rand_marg(s,2)-delta_del;
            rand_marg(t,2) = rand_marg(t,2)+delta_del;

            rand_margs(c,s,1) = rand_margs(c,s,1)-delta_amp;
            rand_margs(c,s,2) = rand_margs(c,s,2)-delta_del;
            rand_margs(c,t,1) = rand_margs(c,t,1)+delta_amp;
            rand_margs(c,t,2) = rand_margs(c,t,2)+delta_del;
            temp = idx_mat(c,s);
            idx_mat(c,s) = idx_mat(c,t);
            idx_mat(c,t) = temp;
            a = a+1; % count backstep
        end
%        end
    end % sub-iterations at same temperature & lineage
    
    % time to adjust temperature?
    if mod(l ,ceil(g/opts.ramp_steps)) == 0
        % adjust temperature, save copy of random marginals
        Temp_resets = Temp_resets + 1;
        Temp = Temp .* step;
    end

    % time to log statistics?
    if mod(l ,opts.stat_interval) == 0
        data.stats(pp,:) = sum(abs((rand_marg-marg)./max(marg,2500))); %!
        data.temp(pp,:) = Temp;
        data.ratios(pp) = k;
        k = 0;
        data.backsteps(pp) = a;
        a = 0;
        pp = pp+1;
    end
end

stat_final = sum(abs((rand_marg-marg)./max(marg,2500))); %! switch from L2 to L1 readout 
%!stat_final = sum(((rand_marg-marg)./max(marg,2500)).^2);


end

