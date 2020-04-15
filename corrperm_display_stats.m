function h = corrperm_display_stats(stats)

posx = [.1,.95];
posy = [.1,.30,.35,.95];

h = figure;
h1 = subplot('Position',[posx(1),posy(3),posx(2)-posx(1),posy(4)-posy(3)]);
% plot error
hold on
for i=1:length(stats);
    N = size(stats{i}.stats,1);
    plot(1:N,stats{i}.stats(:,1),'r',1:N,stats{i}.stats(:,2),'b');
end
ylabel('relative square \Delta')
set(gca,'YScale','log')
set(gca,'TickDir','out')
set(gca,'XTickLabel',[])

% plot temperature
h2 = subplot('Position',[posx(1),posy(1),posx(2)-posx(1),posy(2)-posy(1)]);
hold on
for i=1:length(stats);
    N = size(stats{i}.temp,1);
    plot(1:N,stats{i}.temp(:,1),'r',1:N,stats{i}.temp(:,2),'b');
end
set(gca,'YScale','log')
ylabel('1/kT')
xlabel('time')
set(gca,'TickDir','out')

linkaxes([h1,h2],'x');