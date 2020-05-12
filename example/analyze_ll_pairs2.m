%% analyze low-level pair correlations - fragment

options = struct;
options.sig_thresh = 0.95;
options.power_thresh = 0.1;
options.tail = 'both';
options.split_eq = true;
options.ext = 'll_pair';
options.lineage_out = unique({D.sis.disease});
options.perm_file_mask = 'idx_cell.chunk.*.mat';
options.pcount = 0;

corrperm_analyze_pairs2(ref_dir,perm_dir,results_dir,options);


