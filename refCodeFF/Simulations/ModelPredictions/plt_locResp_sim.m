function plt_locResp_sim(p_r, s, model, s_V, ttl, bool_save, bool_plt_sim_r, counts)
%by defualt, don't save the plot and don't plot the simulated r (only plot
%p(r_A|s_A, s_V) and p(r_V|s_A, s_V)
if nargin < 7; bool_plt_sim_r = 0; counts = [];end
if nargin < 6; bool_save = 0;end
%define plotting info
global lenS lenM
modality_sign = [1,-1];
cb            = [65,105,225]./255; cr = [0.85, 0.23, 0.28];
cMAP          = {cb, cr};
x_bds         = [min(s), max(s)];
y_bds         = [-0.2, 0.2];
x_ticks       = [min(s), s_V(1:2:end), max(s)];
y_ticks       = -0.125:0.125:0.125;
lw            = 2;

figure
for i = 1:lenS
    subplot(1,lenS, i)
    addBackground(x_bds, y_bds, x_ticks, y_ticks);
    for n = 1:lenM
        %fill the area under the distributions
        patch([model.bins_r, fliplr(model.bins_r)], ...
            [modality_sign(n).*squeeze(p_r(i,n,:))',...
            zeros(1,model.numBins_r)],cMAP{n},'FaceAlpha',0.2,...
            'EdgeColor',cMAP{n},'EdgeAlpha',0.9,'LineWidth',2); hold on;
        %plot the MAP distributions
        if bool_plt_sim_r
            plot(model.bins_r, modality_sign(n).*squeeze(counts(i,n,:))',...
                'Color', cMAP{n},'lineWidth',1.5); hold on;
        end
    end
    %plot the mean of biased estimates
    h(1) = plot([s_V(1), s_V(1)],[0, y_bds(end)],'Color',cb,'lineWidth',lw,...
            'lineStyle','--'); hold on;
    h(2) = plot([s_V(i), s_V(i)],[y_bds(1),0],'Color',cr,'lineWidth',lw,...
            'lineStyle','--'); hold on;
    %the horizontal axis that divides modalities
    plot(s, zeros(1,model.numBins_r),'Color', [0,0,0],'lineWidth',lw,...
        'lineStyle','-'); hold off;box off;
    %tick labels and tick marks
    if i == 1
        ylabel('p(R=r|s_{AV,A},s_{AV,V})'); 
        yticks(y_ticks);yticklabels({'0.125','0','0.125'}); 
    else; yticks([]);
    end
    if i == ceil(lenS/2); xlabel('Stimulus location (dva)');end
    xticks(s_V(i));
    %legend
    if i == 1
        legend([h(1), h(2)], {'Auditory','Visual'}, 'Location',...
            'northeast','FontSize',20); legend boxoff;
    end
    xlim(x_bds); ylim(y_bds); set(gca,'FontSize',25);
end
%add title        
sgtitle(ttl, 'FontSize',20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.4]);
set(gcf,'PaperUnits','centimeters','PaperSize',[65 17]);
if bool_save; saveas(gcf, ['LocResp_sim_',ttl], 'pdf'); end

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

