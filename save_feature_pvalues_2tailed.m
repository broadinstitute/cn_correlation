function save_feature_pvalues_2tailed(regs,features,list,filename,p_to_q,sig_thresh)
% SAVE_PVALUES_2TAILED - peak-vs-peak results, pval corrected for 2 tails
%
%   save_feature_pvalues(REGS,FEATURES,LIST,FILENAME,P_TO_Q,SIG_THRESH)
%
% REGS - 1 X 2 cell array containing amp/del peaks regions w/'name' field added
% FEATURES - label for each sample feature
% LIST - Npairs X 3 matrix: column 1 is p-value, col 2: 1=correlated/2=anti ,
%        col 3 is a peak index, col 4 is a feature index
% FILENAME - name of output file
% P_TO_Q - true to convert p-values to q-values, false to leave as p-values
% SIG_THRESH - uppper limit of p- or q-value to include pair in putput

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

% write output
corr_anti = {'corr';'anti'};
write_filtered_tabcols(filename,[],...
    {'index1',sigindex(siglist(:,3))},...
    {'index2',siglist(:,4)},...
    {'p_valueX2',siglist(:,1)},...
    {'q_value',q},...
    {'scna_label',labels(siglist(:,3))},...
    {'tail',corr_anti(siglist(:,2))},...
    {'feature_label',features(siglist(:,4))} ...
);

function O=grade(val)
[~,O] = sort(val);
