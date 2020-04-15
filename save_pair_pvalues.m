function save_pair_pvalues(regs,list,filename,p_to_q,sig_thresh,power_thresh)
% SAVE_PAIR_PVALUES - output peak-vs-peak results
%
%   save_pair_pvalues(REGS,LIST,FILENAME,P_TO_Q,SIG_THRESH)
%
% REGS - 1 X 2 cell array containing amp/del peaks regions w/'name' field added
% LIST - Npairs X 3 matrix: column 1 is p-value, cols 2 and 3 are peak
%        indices corresponding to REGS (amps followed by dels). An optional
%        fourth column may hold the p-value for a hypothesis power filter.
% FILENAME - name of output file
% P_TO_Q - true to convert p-values to q-values, false to leave as p-values
% SIG_THRESH - uppper limit of p- or q-value to include pair in putput
% POWER_THRESH - upper limit on power p-value (optional filter parameter)

if ~exist('power_thresh','var') || isempty(power_thresh)
    power_thresh = 1;
end

% calculate index by significance for each peak (like NG paper)
sigindex = [grade(grade([regs{1}.resid_qv])),grade(grade([regs{2}.resid_qv]))+length(regs{1})];

% optionally filter list to keep only powered hypotheses
if power_thresh < 1 && size(list,2) > 3
    pkeepers = list(:,4) <= power_thresh;
    list = list(pkeepers,:);
end

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

% if max power provided as input, pass it to output.
if size(list,2) > 3
    power_column = {{'max.power',siglist(:,4)}};
else
    power_column = {};
end

% write output
output_cols = { ...
    {'index1',sigindex(siglist(:,2))},...
    {'index2',sigindex(siglist(:,3))},...
    power_column{:},...
    {'p_value',siglist(:,1)},...
    {'q_value',q},...
    {'event_1_label',labels(siglist(:,2))},...
    {'event_2_label',labels(siglist(:,3))} ...
};
write_filtered_tabcols(filename,[],output_cols{:});


function O=grade(val)
[~,O] = sort(val);
