function plotVE(VE,modality,SEM,spatialD)
cb = [65,105,225]./255; cr = [0.85, 0.23, 0.28];
cMap = [cb;cr];
x_bds = [spatialD(1)-2, spatialD(end)+2]; y_bds = [-2,9];
x_ticks = spatialD; y_ticks = 0:3:10;

figure
addBackground(x_bds, y_bds, [x_bds(1), x_ticks', x_bds(end)], y_ticks)
for i = 1:length(modality)
    h(i) = errorbar(spatialD, VE(i,:), SEM(i,:),'-o',...
        'MarkerSize',20,'Color',cMap(i,:),...
        'MarkerFaceColor',cMap(i,:),'MarkerEdgeColor',...
        cMap(i,:),'lineWidth',3); hold on
    hold on;
end
box off; legend([h(1), h(2)],{'Auditory', 'Visual'},'Location','northwest'); 
legend boxoff; xticks(x_ticks); yticks(y_ticks);
xlim(x_bds); ylim(y_bds); xlabel('Absolute spatial discrepancy (dva)');
ylabel('Averaged shifts (dva)'); set(gca,'FontSize',25);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.32, 0.55]);
set(gcf,'PaperUnits','centimeters','PaperSize',[25 25]);
saveas(gcf, 'STATS_VE', 'pdf'); 

%Add background (similar to Seaborn)
function addBackground(x_bd, y_bd, x_ticks, y_ticks)
patch([x_bd(1), x_bd(2), x_bd(2), x_bd(1)], [y_bd(1), y_bd(1), y_bd(2), y_bd(2)],...
    [234,234,242]./255, 'EdgeAlpha',0); hold on
if length(xticks) >=3 
    for i = 2:length(x_ticks)-1
        plot([x_ticks(i), x_ticks(i)], [y_bd(1), y_bd(2)], 'Color',...
            [1,1,1],'lineWidth', 0.5); hold on;
    end
end
if length(yticks) >= 3
    for j = 2:length(y_ticks)-1
        plot([x_bd(1), x_bd(2)],[y_ticks(j), y_ticks(j)], 'Color',...
            [1,1,1],'lineWidth', 0.5); hold on;
    end
end
