function plot_regLines(a_A, b_A, lb_a_A, ub_a_A, lb_b_A, ub_b_A, mean_a_A,...
    mean_b_A, SE_a_A, SE_b_A, bool_save)
x_bds    = [0, 1.5]; y_bds = [-3, 12];
x_ticks  = 0:0.5:1.5; y_ticks = 0:3:9;

figure
addBackground(x_bds, y_bds, x_ticks, [y_bds(1), y_ticks, y_bds(end)])
plot(x_bds,[0,0],'k--','lineWidth',1.5); hold on;
plot([1,1],y_bds,'k--','lineWidth',1.5); hold on;
h1=errorbar(a_A, b_A, b_A - lb_b_A, ub_b_A - b_A,a_A - lb_a_A, ub_a_A - a_A,...
    'o','MarkerSize',10,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',...
    0.6.*ones(1,3),'Color',0.6.*ones(1,3),'lineWidth',1.5,'CapSize',10); hold on;
h2 = errorbar(mean_a_A, mean_b_A, SE_b_A, SE_b_A,SE_a_A, SE_a_A, 'o',...
    'MarkerSize',10,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[255,255,102]./255,...
    'Color',[0,0,0],'lineWidth',2,'CapSize',10); hold off; box off;
xlim(x_bds); ylim(y_bds); xticks(x_ticks); yticks(y_ticks);
xlabel(sprintf('The slope of the linear regression\nfor PSEs on visual stimulus location'));
ylabel(sprintf('The intercept of the linear regression\nfor PSEs on visual stimulus location'));
legend([h1,h2],{'Individual data', 'Group data'},'FontSize',15,'Location','northwest'); 
legend boxoff; set(gca,'FontSize',18);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.3, 0.55]);
set(gcf,'PaperUnits','centimeters','PaperSize',[20 25]);
if bool_save; saveas(gcf, 'RelativeBiases_allSubjs', 'pdf'); end

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