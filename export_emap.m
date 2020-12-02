function export_emap(E,evif_fname,ecall_fname)
%EXPORT_EMAP write event map as tab-delimited text files
%
%    export_emap(EMAP,EVIF_FNAME,ECALL_FNAME)
%
% EMAP is an event mapping structure used to analyze permutations
% EVIF_FNAME is the name of the tab-delimited event information file to write
% ECALL_FNAME is the name of the tab-delimited event call file to write
%
    %% write event map text file
    tb = char(9); % tab character delimiter
    fad = {'f','a','d'}; % map integer to event type
    O12 = {'0','1','2'}; % map logical matrix to char

    % write event info file
    other_fields = setdiff(fieldnames(E.event),{'name','type','chrn'});
    other_cols = {};
    for f = 1:length(other_fields)
        event_data = E.event.(other_fields{f});
        if isnumeric(event_data) && any(event_data-round(event_data))
            other_cols{f} = {other_fields{f},E.event.(other_fields{f}),'%0.16g'};
        else
            other_cols{f} = {other_fields{f},E.event.(other_fields{f})};
        end
            
    end

    write_filtered_tabcols(evif_fname,[],{'event',E.event.name},...
                           {'type',fad(E.event.type+1)},...
                           {'chrn',E.event.chrn},other_cols{:});

    % write event call file 
    fid = fopen(ecall_fname,'w');
    % construct 'lineage' output class column
    lincol = repmat({''},size(E.sample.id));
    for i = 1:length(E.pcx)
        lincol(E.pcx{i}) = repmat(cellstr(E.pcname{i}),size(E.pcx{i}));
    end
    % write header
    fprintf(fid,['sample',tb,'lineage',tb,strjoin(E.event.name',tb),'\n']);
    % write body
    for i = 1:size(E.dat,2)
        fprintf(fid,[E.sample.id{i},tb,lincol{i},tb,strjoin(O12(E.dat(:,i)'+1),tb),'\n']);
    end
    fclose(fid);
end
