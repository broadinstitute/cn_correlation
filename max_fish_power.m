function [pmin_corr,pmin_anti] = max_fish_power(Nsamples,cd,bd)
%MAX_FISH_POWER minimum Fisher exact p-values possible, given marginals
%
%   [PMIN_CORR,PMIN_ANTI] = max_fish_power(NSAMPLES,E1,E2)
%

%        key  ~E2   E2
%         ~E1  a    b
%          E1  c    d
%
        % correlation (maximum a,d)
        d_mc = min(cd,bd);              % maximum coincidence
        b_mc = bd - d_mc;
        c_mc = cd - d_mc;
        a_mc = Nsamples - b_mc - c_mc - d_mc;
        pmin_corr = fisher_exact_test(a_mc,b_mc,c_mc,d_mc);

        % anti-correlation (maximum b,c)
        ab = Nsamples - cd;             % reformulate constraint
        b_ma = min(ab,bd);
        a_ma = ab - b_ma;
        d_ma = bd - b_ma;
        c_ma = cd - d_ma;
        pmin_anti = fisher_exact_test(a_ma,b_ma,c_ma,d_ma);
