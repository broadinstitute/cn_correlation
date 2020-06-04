function export_emap(E,fname)
%EXPORT_EMAP write event map to a text file
%
%    export_emap(EMAP,FNAME)
%
% EMAP is an event mapping structure used to analyze permutations
% FNAME is the path to the .txt file output
%
    %% write event map text file
    tb = char(9); % tab character delimiter
    fad = {'f','a','d'}; % map integer to event type
    O12 = {'0','1','2'}; % map logical matrix to char
    fid = fopen(fname,'w');
    fprintf(fid,['name',tb,strjoin(E.event.name',tb),'\n']);
    fprintf(fid,['type',tb,strjoin(fad(E.event.type+1),tb),'\n']);
    fprintf(fid,['chrn',tb,strjoin(cellstr(num2str(E.event.chrn))',tb),'\n']);
    for i = 1:size(E.dat,2)
        fprintf(fid,[E.sdesc{i},tb,strjoin(O12(E.dat(:,i)'+1),tb),'\n']);
    end
    fclose(fid);
end
