function output_pair_p(events,list,filename,p_to_q,sig_thresh,power_thresh)
% OUTPUT_PAIR_P - output peak-vs-peak results
%
%   output_pair_p(EVENTS,LIST,FILENAME,P_TO_Q,SIG_THRESH)
%
% EVENTS - struct of per-event arrays, relevant fields are 'name' and 'resid_qv'
% LIST - Npairs X 3 matrix: column 1 is p-value, cols 2 and 3 are event
%        indices corresponding to EVENT. An optional fourth column
%        may hold the p-value for a hypothesis power filter.
% FILENAME - name of output file
% P_TO_Q - true to convert p-values to q-values, false to leave as p-values
% SIG_THRESH - uppper limit of p- or q-value to include pair in putput
% POWER_THRESH - upper limit on power p-value (optional filter parameter)

    if ~exist('power_thresh','var') || isempty(power_thresh)
        power_thresh = 1;
    end

    % calculate index by significance for each peak within event type (GISTIC peaks output order)

    %! amps = events.type == 1;
    %! dels = events.type == 2;
    %! sigindex = [grade(grade(events.resid_qv(amps)')), grade(grade(events.resid_qv(dels)')) + sum(amps)];

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

    %!labels = [strcat('amp_',{regs{1}.name}),strcat('del_',{regs{2}.name})];
    labels = events.name;

    % if max power provided as input, pass it to output.
    if size(list,2) > 3
        power_column = {{'max.power',siglist(:,4)}};
    else
        power_column = {};
    end

    % write output
    output_cols = { ...
        {'event_1_label',labels(siglist(:,2))},...
        {'event_2_label',labels(siglist(:,3))},...
        power_column{:},...
        {'p_value',siglist(:,1)},...
        {'q_value',q}...
                  };
    write_filtered_tabcols(filename,[],output_cols{:});
end

%!function O=grade(val)
%!    [~,O] = sort(val);
%!end

