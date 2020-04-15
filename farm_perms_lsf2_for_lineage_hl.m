% compile instructions
% mbuild -setup % choose long answer
% mcc -v -m -w enable annealing_permutations_lsf2_hl.m
%!use -q LSF
N = 100; % number of jobs to run
M = 50; % number of iterations
command = './permwrap_lsf2_hl.sh';
%!command = './annealing_permutations_lsf2';

input_dir = [pwd,'/input_hl_0.2_sub/'];
output_dir = [pwd,'/output_hl_0.2_sub/'];
chrm_metro_restart_lsf2_lineage_hl(D,regs,input_dir)

cmd_template = ['bsub ',...
          '-R "rusage[mem=4]" ',... 
          '-mig 5 ',...
          '-R "select[cpuf>100]" ',... 
          '-Q "EXCLUDE(127)" ',...
          '-q hour ',...
          '-P cancerfolk ',...
          '-o ',output_dir,'%s.out.txt -e ',output_dir,'%s.err.txt ',...
          '-r ',command,' ',input_dir,' ',output_dir ' %s %d']; 
      
      
% need input directory
if ~exist(input_dir,'dir')
    error('input directory doen''t exist!')
end
% make output directory if need be
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end
% submit to the farm
for j=1:N
    runid = sprintf('cycle%03d',j);
    unix_str = sprintf(cmd_template,runid,runid,runid,M);
    disp(unix_str);
   [r1,r2]=unix(unix_str); 
   disp([strtrim(r2) ' Exit code: ' num2str(r1)])
end