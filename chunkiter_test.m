% concrete example of calling chunkiter
set_verbose_level(50)
addpath ~/git/broadinstitute/cn_correlation
addpath ~/git/broadinstitute/cn_correlation/snputil

perm_dir = '/xchip/beroukhimlab/gistic/corrperm/workspace/ll_work/ll_permout';
emap_file = '/xchip/beroukhimlab/gistic/corrperm/example/ll_work/ngll.emap.mat';
options = struct('analyze_lineages',true);

load(emap_file);
E.dat = full(E.dat);
[ot,mlt] = chunkiter(E,perm_dir,options);
