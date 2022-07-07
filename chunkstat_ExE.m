% chunk-processor implementation for standard co-ocurrence p-value processing
classdef chunkstat_ExE < chunkstat
    properties
        le_counts = 0;
        ge_counts = 0;
        eq_counts = 0;
        obs_tot = 0;
    end

    methods
        % constructor called before chunk processing
        % S = chunkstat(emap)
        % emap is an Nevents x Nsamples logical array of observed events in each sample
        function S = chunkstat_ExE(emap)
            emap = double(emap);
            S@chunkstat(emap)
            % ccount observed co-occurences
            S.obs_tot = emap * emap';

            %![Nevents,Nsamples] = size(E.dat);
            
            % initialize E x E count
            S.le_counts = zeros(S.Nevents,S.Nevents);
            S.ge_counts = zeros(S.Nevents,S.Nevents);
            S.eq_counts = zeros(S.Nevents,S.Nevents);
        end
        
        %% called for every permutation
        function S = perm(S,emap)
            % call superclass
            S = perm@chunkstat(S,emap);

            % count co-occurrences
            cooc = emap * emap';

            % compare with observed, count equality and each direction of inequality
            S.le_counts = S.le_counts + (cooc <= S.obs_tot);
            S.ge_counts = S.ge_counts + (cooc >= S.obs_tot);
            S.eq_counts = S.eq_counts + (cooc == S.obs_tot);
        end

        % called after chunk processing complete
        function table = results(S,pairs_idx,options)
            % allocate storage for lineage-specific p-values
            Npairs = size(pairs_idx,1);   % peak indices
            p_corr = zeros(Npairs,1);
            p_anti = zeros(Npairs,1);

            % calculate p-values

            verbose('calculating overall p-values for co-occurrences',20);
            for s = 1:Npairs
                i = pairs_idx(s,1);
                j = pairs_idx(s,2);
                
                % p value calculations
                pscnt = options.pcount;
                N = S.Nperms + pscnt;
                if (options.split_eq)
                    p_corr(s) = (S.ge_counts(i,j) - S.eq_counts(i,j)/2 + pscnt) / N;
                    p_anti(s) = (S.le_counts(i,j) - S.eq_counts(i,j)/2 + pscnt) / N;
                else
                    p_corr(s) = (S.ge_counts(i,j) + pscnt) / N;
                    p_anti(s) = (S.le_counts(i,j) + pscnt) / N;
                end

            end
            table = struct('p_corr',num2cell(p_corr),...
                           'p_anti',num2cell(p_anti) );
        
        end % function
        
    end % methods
end
