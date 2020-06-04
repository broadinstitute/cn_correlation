% EDIT_ME script to run a few permutations on local CPU and display a plot of runs for schedule tuning
% uses global variables H, perm_opts
tune_trials = 5;
[~,stat_finals,stats] = ampdel_tempering(H,tune_trials,perm_opts);
corrperm_display_stats(stats);
figure; % save plot of trial run schedule
corrperm_display_stats(stats);
save_current_figure(fullfile(pwd,'ll.tuning_stats'),{'png','pdf'});
% display mean, spread of error
mean_err = mean(vertcat(stat_finals{:}))
std_err = std(vertcat(stat_finals{:}))
