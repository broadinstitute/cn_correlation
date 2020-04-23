function structable = read_R_table(fname,dlm,nantext,options)
% READ_R_TABLE read an R-style table as a Matlab struct
%
%  STRUCTABLE = read_R_table(FNAME,DLM,NANTEXT,OPTIONS)
%
% Read an R-style delimited text file FNAME containg a data table with a
% header into a MATLAB struct, with the field names formed from the column
% headers. DLM specifies the delimiter, default tab=char(9). NANTEXT 
% specifies the token used for missing numeric data, default 'NA'. Missing
% numeric data will be signified by NaN in the struct array data.

% default tab-delimited table
if ~exist('dlm','var') || isempty(dlm)
    dlm = char(9);
end
% 
if ~exist('nantext','var')
    nantext = 'NA';
end

if ~exist('options','var')
    % structure of options can be at any position of legacy inputs
    if nargin==2 && isa(dlm,'struct')
        options = dlm;
        dlm = char(9);
    elseif nargin==3 && isa(nantext,'struct')
        options = nantext;
        nantext = 'NA';
    else
        options = struct();
    end
end

options = impose_default_value(options,'colClasses','');
%!!! to really be like R, should implement colClasses as vector
options = impose_default_value(options,'sep',dlm);
options = impose_default_value(options,'nantext',nantext);
nantext = ['^',options.nantext,'$'];

% read header for number of columns
fid = fopen(fname);
skipline = true;
while skipline
    hdrline = fgets(fid);
    skipline = hdrline(1) == '#';
end
hdr = regexp(hdrline,options.sep,'split');
Ncols = length(hdr);

% make textscan format string
if hdrline(end) == char(10)
    lf_char = '\n';%char(10);
    pcunix_text = true;
else
    lf_char = '\r';%char(13);
    pcunix_text = false;
end
fmstr = [repmat('%s ',1,Ncols-1),'%s',lf_char];

% read in file as cell array of columns
try
    if pcunix_text
        cols = textscan(fid,fmstr,'ReturnOnError',0,'Delimiter',char(9),'WhiteSpace','\r');
    else
        cols = textscan(fid,fmstr,'ReturnOnError',0,'Delimiter',char(9));
    end
    fclose(fid);
catch me
    fclose(fid);
    throw(MException('snp:error',...
                    strcat('Error scanning segment file ''%s'':\n',me.message),...
                    fname));
end

% pad any short columns (while swearing at excel)
colens = cellfun(@length,cols);
if any(diff(colens))
    maxlen = max(colens);
    for c = 1:length(colens);
        if colens(c) < maxlen
            cols{c} = [cols{c};repmat({''},maxlen-colens(c),1)];
        end
    end
end

cpp = horzcat(cols{:});

%{
cpp = read_dlm_file(fname,dlm);
comment_lines = cellfun(@(x) ~isempty(regexp(x{1},'^#','start','once')),cpp);
if any(comment_lines)
    cpp{comment_lines} = [];
end
cpp = vertcat(cpp{:});
%}

if isempty(cpp)
    structable = [];
else

    % convert header into legitimate structure field name
    hdr = regexprep(strtrim(hdr),'\s+','_'); %!
%!  hdr = regexprep(strtrim(cpp(1,:)),'\s+','_');
    hdr = regexprep(hdr,'[^A-Za-z0-9_]+','');
    hdr = regexprep(hdr,'^([0-9_]+)','X$1');
    
    for i = 1:length(hdr)-1
        dups = find(strcmp(hdr{i},hdr(i+1:end)));
        for j = 1:length(dups)
            dup = dups(j);
            n = 1;
            while any(strcmp([hdr{i+dup},'_',num2str(n)],hdr))
                n = n + 1;
            end
            hdr{i+dup} = [hdr{i+dup},'_',num2str(n)];
        end
    end

    ncols = length(hdr);
    nrows = size(cpp,1);%!-1;
    structcell = cell(2,ncols);
    structcell(1,:) = hdr;
    for col = 1:ncols
        if strcmp(options.colClasses,'char')
            structcell{2,col} = cpp(:,col);
        else
            % attempt to convert to number
            trynum = str2num(char(regexprep(cpp(:,col),nantext,'NaN')));
            if length(trynum) == nrows
                structcell{2,col} = num2cell(trynum);
            else
                % cannot be made numeric
                structcell{2,col} = cpp(:,col);
            end
        end
    end
    structable = struct(structcell{:});
end
