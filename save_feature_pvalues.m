function save_feature_pvalues(regs,features,list,filename,p_to_q,sig_thresh,power_thresh)
% SAVE_PAIR_PVALUES - output peak-vs-feature results
%
%   save_feature_pvalues(REGS,FEATURES,LIST,FILENAME,P_TO_Q,SIG_THRESH)
%
% REGS - 1 X 2 cell array containing amp/del peaks regions w/'name' field added
% FEATURES - label for each sample feature
% LIST - Npairs X 3 matrix: column 1 is p-value, col 2 is a peak index,
%        col 3 is a feature index. A fourth column may be added to provide
%        information about the experimental power for the hypothesis, a
%        fifth column may be provided to indicate which values are
%        correlations (nonzero) and which are anti-correlations (zero)
% FILENAME - name of output file
% P_TO_Q - true to convert p-values to q-values, false to leave as p-values
% SIG_THRESH - uppper limit of p- or q-value to include pair in putput
% POWER_THRESH - upper limit on power p-value (optional filter parameter)

if ~exist('power_thresh','var') || isempty(power_thresh)
    power_thresh = 1;
end

% optionally filter list to keep only powered hypotheses
if power_thresh < 1 && size(list,2) > 3
    pkeepers = list(:,4) <= power_thresh;
    list = list(pkeepers,:);
end

% calculate index by significance for each peak (like NG paper)
sigindex = [grade(grade([regs{1}.resid_qv])),grade(grade([regs{2}.resid_qv]))+length(regs{1})];

% filter list by p- or q-value
if p_to_q
    q = calc_fdr_value(list(:,1));
    keepers = q<=sig_thresh;
else
    keepers = list(:,1)<=sig_thresh;
end
siglist = list(keepers,:);
q = q(keepers);
labels = [strcat('amp_',{regs{1}.name}),strcat('del_',{regs{2}.name})];

%% write output

% if max power provided as input, pass it to output.
if size(list,2) > 3
    power_column = {{'max.power',siglist(:,4)}};
else
    power_column = {};
end

% if correlation/anticorrelation provided as input, indicate in output
if size(list,2) > 4
    corr_anti = {'corr','anti'};
    ac_column = {{'tail',corr_anti(1 + ~siglist(:,5))}};
else
    ac_column = {};
end

output_cols = { ...
    {'scna_index',sigindex(siglist(:,2))},...
    {'feature_index',siglist(:,3)},...
    power_column{:},...
    ac_column{:},...
    {'p_value',siglist(:,1)},...
    {'q_value',q},...
    {'scna_label',labels(siglist(:,2))},...
    {'feature_label',features(siglist(:,3))} ...
};

% write tab-delimited columns to file
write_filtered_tabcols(filename,[],output_cols{:});

function O=grade(val)
[~,O] = sort(val);
