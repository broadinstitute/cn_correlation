% EDIT_ME script to submit permutations to multiprocessor environment
% uses global variables ref_dir, perm_dir, perm_opts
submit_perms(H,[perm_opts.algorithm,'_module'],'uger.submit',ref_dir,perm_dir,100,50,perm_opts);
