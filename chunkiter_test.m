% concrete example of calling chunkiter
set_verbose_level(50)
addpath ~/git/broadinstitute/cn_correlation
addpath ~/git/broadinstitute/cn_correlation/snputil

output_dir = '/xchip/beroukhimlab/gistic/corrperm/workspace/gew_2017';
perm_dir = '/xchip/beroukhimlab/gistic/corrperm/workspace/gew_2017/permout';
emap_file = '/xchip/beroukhimlab/gistic/corrperm/workspace/gew_2017/gew_gg.score_D_scnas.emap.mat';
sif_file = '/xchip/beroukhimlab/gistic/corrperm/workspace/gew_2017/gew.160229.sif.txt';
options = struct('analyze_lineages',true);

E = load_emap(emap_file);

% add is_west feature in a tortured fashion
SI = read_R_table(sif_file);

[~,sx,ex] = match_string_sets_hash({SI.sample},{E.sample.id});
is_west = num2cell(strcmp({SI(sx).implied_race},'west'));
[E.sample(ex).is_west] = deal(is_west{:});
E.dat = [E.dat;[is_west{:}]];
E.event.name = [E.event.name;'is_west'];
E.event.chrn = [E.event.chrn;0];
E.event.type = [E.event.type;0];
E.event.resid_qv = [E.event.resid_qv;NaN];

results = chunkiter(E,perm_dir,options);

names = fields(results);
for f = 1:length(names)
    savestruct(results.(names{f}),fullfile(output_dir,['results.',names{f},'.txt']));
end
