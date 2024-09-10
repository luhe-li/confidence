function plot_PSE(D, Expt, Plt)
vLoc = [-12,-4,4,12];
cMAP = [255,178,178;255,127,127;229,0,0;153,0,0]./255;
x_bds = [vLoc(1)-5, vLoc(end)+5]; y_bds = x_bds;
x_ticks = vLoc; y_ticks = round(D.polyval,1);

figure
addBackground(x_bds, y_bds, [x_bds(1), x_ticks, x_bds(end)], [y_bds(1), y_ticks, y_bds(end)]);
%plot the 4 PSE's with error bars
for i = 1:Expt.testLocations
    errorbar(vLoc(i), D.estimatedP(2*i),D.estimatedP(2*i)-...
        D.PSE_lb(i),D.PSE_ub(i)-D.estimatedP(2*i),'s',...
        'MarkerSize',20,'MarkerEdgeColor',cMAP(i,:),'MarkerFaceColor',...
        cMAP(i,:),'Color',cMAP(i,:),'LineWidth',3,'CapSize',10);
    hold on
end
%plot the fitted straightline
h1 = plot(vLoc,D.polyval,'k-','LineWidth',2); hold on
%plot the 4 test locations
h2 = plot(vLoc,vLoc,'k--','LineWidth',2); hold on; 
text(x_bds(1)+0.25,y_bds(end)-1,Plt.subjI,'fontSize',15); hold off; box off
xlabel('Location of visual standard stimulus (dva)'); xlim(x_bds);
ylabel('Point of subjective equality (dva)');  ylim(y_bds);
legend([h1, h2],{'Linear regression','Unity line'},'FontSize',15,'Location','southeast'); legend boxoff
set(gca,'FontSize',18, 'YTick', y_ticks, 'XTick', x_ticks); set(gca,'FontSize',18); 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.21, 0.4]);
if Plt.bool_save == 1
    set(gcf,'PaperUnits','centimeters','PaperSize',[15 15]);
    saveas(gcf, ['PSE_sub', num2str(Expt.subjID)], 'pdf'); 
end

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