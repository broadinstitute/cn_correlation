function E = import_emap(fname)
%IMPORT_EMAP read an event map from a text file
%
%     EMAP = import_emap(FNAME)
%
% EMAP is an event mapping structure used to analyze permutations
% FNAME is the path to the .txt file output
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

dlm = char(9); % default tab-delimited table
nantext = 'NA'; %!!! disallow nans, eliminate
nantext = ['^',nantext,'$'];

%% read header for number of columns
fid = fopen(fname);

% scan over empty or R-style # comment lines
skipline = true;
while skipline
    hdrline = fgets(fid);
    skipline = isempty(hdrline) || hdrline(1) == '#';
    %!    fprintf(hdrline); %!!!debug
end

hdr = regexp(hdrline,dlm,'split');
Ncols = length(hdr);
if Ncols == 0
    error('event.map error: empty file');    
end
try_name = strtrim(hdr(1));
if ~strcmp(try_name,'name')
    error('event.map error: first row/column of event map must be ''name''.');
end
Nevents = Ncols-1;
event = struct;
event.name = strtrim(hdr(2:Ncols)');

% process header line attributes
no_chrn_yet = true;
anames = cell(0); % header/attribute names so far
while no_chrn_yet
    anames = [anames,{try_name}];
    hdrline = fgets(fid);
    hdr = regexp(hdrline,dlm,'split');
    if (length(hdr) ~= Ncols)
        error('event.map error: uneven column count');
    end
    try_name = strtrim(hdr(1));
    try_vals = strtrim(hdr(2:Ncols)');
    if strcmp(try_name,'chrn')
        event.chrn = chromosome2num(try_vals);
        no_chrn_yet = false;
    else 
        if strcmp(try_name,'type')
            try_nums = str2num(char(regexprep(regexprep(regexprep(try_vals,'f','0'),'a','1'),'d','2')));
            if length(try_nums) == Nevents
                event.type = try_nums;
            else
                error('event.map error: bad event ''type'' values.');
            end
            
        else % unrecognized attribute
             % clean up name to be valid field name
            try_name = regexprep(strtrim(try_name),'\s+','_'); %!
            try_name = regexprep(try_name,'[^A-Za-z0-9_]+','');
            try_name = regexprep(try_name,'^([0-9_]+)','X$1');
            % name de-duplication
            if any(strcmp(try_name,anames))
                n = 2;
                try_n = try_name;
                while any(strcmp(try_name,anames))
                    try_name = [try_n,'_',num2str(n)];
                    n = n+1;
                end
            end
             % attempt to convert to number
            try_nums = str2num(char(regexprep(tryvals,nantext,'NaN')));
            if length(try_nums) == Nevents
                % save as number
                event.(try_name) = try_nums;
            else
                % keep as string
                event.try_num = try_vals;
            end
            %    event(hdr(1)) = hdr(2:Ncols)
        end
    end
end

%% read in rest of file as cell array of columns

% make textscan format string
if hdrline(end) == char(10)
    lf_char = '\n';%char(10);
    pcunix_text = true;
else
    lf_char = '\r';%char(13);
    pcunix_text = false;
end

% the first column is sample label, the rest numbers
fmstr = ['%s',repmat(' %d',1,Nevents),lf_char];
try
    if pcunix_text
        cols = textscan(fid,fmstr,'ReturnOnError',0,'Delimiter',dlm,'WhiteSpace','\r');
    else
        cols = textscan(fid,fmstr,'ReturnOnError',0,'Delimiter',dlm);
    end
    fclose(fid);
catch me
    fclose(fid);
    error(['event map error: textscan error ''',me.message,'''']);
end

% test for short columns
if any(diff(cellfun(@length,cols)))
    error('event map error: missing event data');
end

E = struct;
E.sdesc = cols{1};
E.event = event;

Nsamples = length(E.sdesc);
E.dat = false(Nevents,Nsamples);
for i = 1:Nevents
    E.dat(i,:) = 0 ~= cols{i+1};
end

end % function
