function [D,samples,reorder] = create_D_sample_bins(D,thresh,sisfield)
%CREATE_D_SAMPLE_BINS - create permutation classes for a D struct
%
%   [D,SAMPLES,REORDER] = create_D_sample_bins(D,THRESH,SISFIELD)
%
% D is a D-struct of copy number. The samples in the returned D will be
% reordered according to permutation classes.
% THRESH is a minimum number of samples allowed in a class,. Classes
% with insufficient samples will be removed from the resulting D.
% SISFIELD designates a field in the D.sis sample information that
% can be used to define classes for reordering.
% The SAMPLES output is a cell array of per-class index vectors into D.
% The REORDER vector is for reording data structures that were aligned with
% the input D to bring them into alignment with the output D.

% (descended from Travis' make_sample_bins, but reorders D in lineage order)

% defaults
if ~exist('thresh','var') || isempty(thresh)
    thresh = 40;
end
if ~exist('sisfield','var') 
    sisfield = [];
end

% no class case
if isempty(sisfield)
    reorder = 1:length(D.sdesc);
    samples = {reorder};
else
    % get unique lineages and per-sample index to them
    if ischar(D.sis(end).(sisfield))
        [utags,~,tagx] = unique({D.sis.(sisfield)},'legacy');
    else
        [utags,~,tagx] = unique([D.sis.(sisfield)],'legacy');
    end
    % see which ones have enough power 
    keeptags = accumarray(tagx',1) > thresh;
    %!cell_type_used = utags(keeptags);

    % find index for lineage order
    [tagx,reorder] = sort(tagx);
    % retain only lineages with power 
    reorder = reorder(ismember(tagx,find(keeptags)));
    D = reorder_D_cols(D,reorder);

    % build permutation classes 
    samples = cell(1,sum(keeptags));
    s = 1;
    for i = 1:length(utags)
       if keeptags(i)
           samples{s} = find(tagx==i);
    %!     samples{s} = find(strcmp({D.sis.(sisfield)},utags{i}));
           s = s+1;
       else
           if ischar(utags)
              badtag = utags{i};
           else
              badtag = num2str(utags(i));
           end
           warning('create_sample_bins:underpowered_lineage',...
                    'removing permutation class %s== ''%s'' - too few samples!',sisfield,badtag);    
       end 
    end
end