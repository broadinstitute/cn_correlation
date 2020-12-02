function E = import_emap(evif_fname,ecall_fname)
%IMPORT_EMAP read an event map from a text file
%
%     EMAP = import_emap(EVIF_FNAME,EVCALL_FNAME)
%
% EMAP is an event mapping structure used to analyze permutations
% EVIF_FNAME is the tab-delimited event information file name to read
% ECALL_FNAME is the tab-delimited event call file name read
%
%
% Read a tab-delimited hash-commented text file containing an array
% of logical sample X event occurrences into an event map. Event attributes
% are defined by multiple header rows and are named in the first column followed
% tyhe attribute values for each event. The first row must define the event 'name',
% using unique identifiers (e.g. 'amp_MCL1'). The last header row processed must define
% 'chrn' (event chromosome). Subsequent to the 'chrn' row, rows are interpreted 
% as sample-labeled (column 1) vectors of occurrence of events,
%  represented as '0' (no event) or '1' (event).
%
% As an input to the corrperm analysis module, the event map must also have the
% event attribute 'type' (1 or 'a' - amplification; 2 or 'd' - deletion; 
% 0 or 'f' - feature).

    dlm = char(9); % tab-delimited tables

    nantext = 'NA'; %!!! disallow nans, eliminate
    nantext = ['^',nantext,'$'];

    %% validate filename arguments
    if ~exist('evif_fname','var') || isempty(evif_fname)
        error('No event information file argument.');
    end
    if ~exist('ecall_fname','var') || isempty(evif_fname)
        error('No event call file argument.');
    end
    if ~exist(evif_fname,'file')
        error(['Event information file ''',evif_fname,''' doesn''t exist.']);
    end
    if ~exist(ecall_fname,'file')
        error(['Event call file ''',ecall_fname,''' doesn''t exist.']);
    end

    E = struct;
    E.event = struct;

    %%
    %% process event information file
    %%
    etab = read_R_table(evif_fname,dlm);
    % event column
    assert(isfield(etab,'event'));
    E.event.name = {etab.event}';
    % type column
    assert(isfield(etab,'type'));
    [ok,idx] = ismember({etab.type},{'f','a','d'});
    assert(all(ok)); % all recognized type
    E.event.type = idx' - 1;
    % chrn column
    assert(isfield(etab,'chrn'));
    if isnumeric( etab(1).chrn );
        E.event.chrn = [etab.chrn]';
    else
        E.event.chrn = chromosome2num({etab.chrn}');
    end

    % field enumeration
    other_fields = setdiff(fieldnames(etab),{'event','type','chrn'});
    for f = 1:length(other_fields)
        if isnumeric(etab(1).(other_fields{f}));
            E.event.(other_fields{f}) = [etab.(other_fields{f})]';
        else
            E.event.(other_fields{f}) = {etab.(other_fields{f})}';
        end
    end

    %%
    %% process event call file
    %%
    fid = fopen(ecall_fname);
    % scan over empty or R-style # comment lines
    skipline = true;
    while skipline
        hdrline = fgets(fid);
        skipline = isempty(hdrline) || hdrline(1) == '#';
        %!    fprintf(hdrline); %!!!debug
    end
    % read header
    hdr = regexp(hdrline,dlm,'split');
    Ncols = length(hdr);
    if Ncols == 0
        error('event.call error: empty file');    
    end
    if Ncols < 2
        error('event.call error: no events')
    end
    % sample column
    try_name = strtrim(hdr(1));
    if ~strcmp(lower(try_name),'sample')
        error('event.call error: first column of event call table must be ''sample''.');
    end
    %{
    % optional lineage column
    try_lineage = strtrim(hdr(2));
    if strcmp(try_lineage,'lineage')
        have_lineage = true;
        Nevents = Ncols-2;
        lin_fmt = ' %s';
        called_events = strtrim(hdr(3:end));
    else
        have_lineage = false;
        Nevents = Ncols-1;
        lin_fmt = '';
        called_events = strtrim(hdr(2:end));
        if ~ismember(try_lineage,event.name)
            error('event.call error: column 2 of event map must be either an event.info name or ''lineage''.');
        end
    end
    %}
    Nevents = Ncols-1;
    called_events = strtrim(hdr(2:end));

    % discard event.info that does not match event.calls
    efields = fieldnames(E.event);
    matching_info = ismember(E.event.name,called_events);
    if ~all(matching_info)
        extra_info = sum(~matching_info);
        warning([num2str(sum(~matching_info)),' events in event.info without event.calls discarded.']);
        for f = 1:length(efields)
            E.event.(efields{f}) = E.event.(efields{f})(matching_info);
        end     
    end

    % add event.calls without event.info as feature events
    [isec,ia,ib] = intersect(called_events,E.event.name,'stable'); %!!! will error on early ML versions
    if length(isec) < length(called_events)
        unk_evx = find(~ismember(called_events,E.event.name))';
        Nnewfeat = length(unk_evx);
        warning([num2str(Nnewfeat),' event.calls present without event.info converted to features.']);
        % create feature event info
        E.event.name = [E.event.name;called_events(unk_evx)'];
        %! E.event.chrn = [E.event.chrn;zeros(Nnewfeat,1)]; % chromosome "zero"
        %! E.event.type = [E.event.chrn;zeros(Nnewfeat,1)]; % feature type
        restfields = setdiff(efields,'name');
        ib = [ib;unk_evx];
        for f = 1:length(restfields)
            if isnumeric(E.event.(restfields{f}))
                E.event.(restfields{f})(ib) = [E.event.(restfields{f});zeros(Nnewfeat,1)];
            else
                E.event.(restfields{f})(ib) = [E.event.(restfields{f});repmat({''},Nnewfeat,1)];
            end
        end
    end

    %% read in rest of event.call file as cell array of columns

    % make textscan format string
    if hdrline(end) == char(10)
        lf_char = '\n';%char(10);
        pcunix_text = true;
    else
        lf_char = '\r';%char(13);
        pcunix_text = false;
    end


    fmstr = ['%s',repmat(' %d',1,Nevents),lf_char];
    %!    fmstr = ['%s',lin_fmt,repmat(' %d',1,Nevents),lf_char];
    try
        if pcunix_text
            cols = textscan(fid,fmstr,'ReturnOnError',0,'Delimiter',dlm,'WhiteSpace','\r');
        else
            cols = textscan(fid,fmstr,'ReturnOnError',0,'Delimiter',dlm);
        end
        fclose(fid);
    catch me
        fclose(fid);
        error(['event.call error: textscan error ''',me.message,'''']);
    end

    % test for short columns
    if any(diff(cellfun(@length,cols)))
        error('event.call error: missing event data');
    end

    % sample identifier column
    E.sample = struct('id',cols{1});
    %! E.sdesc = strtrim(cols{1})
    if length(sample.id) ~= length(unique(sample.id))
        error('duplicate sample identifiers detected in event calls')
    end
    % optional lineage column
    Nsamples = length(E.sample.id);
    %{    
    if have_lineage
        lintext = strtrim(cols{2});
        E.pcname = unique(lintext,'stable')';
        E.pcx = cellfun(@(x) find(ismember(lintext,x)),E.pcname,'UniformOutput',false);
    else
        E.pcx = {};
        E.pcname = {};
    end
    %}
    % data columns
    E.dat = false(Nevents,Nsamples);
    for i = 1:Nevents
        E.dat(i,:) = 0 ~= cols{i+1};
        %! E.dat(i,:) = 0 ~= cols{i+have_lineage+1};
    end
    E.dat = logical(E.dat);

end % function
