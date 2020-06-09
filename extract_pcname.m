function pcnames = extract_pcname(D,pcx,sisfield)
%EXTRACT_PERMCLASS_NAMES find permutation class names IN D
% DOCME!!!
%
    if ~exist('sisfield')
        sisfield = '';
    end
    % look for permclass fields in sample info structure
    if isempty(sisfield)
        if isfield(D,'sis')
            if isfield(D.sis,'disease')
                sisfield = 'disease';
            else
                if isfield(D.sis,'gcmtype')
                    sisfield = 'gcmtype';
                end
            end
        end
    end
    pcnames = cell(size(pcx));
    for i=1:length(pcx)
        if isempty(sisfield)
            pcnames{i} = sprintf('class%d',i);
        else
            j = pcx{i}(1);
            pcnames{i} = D.sis(j).(sisfield);
        end
    end
end

