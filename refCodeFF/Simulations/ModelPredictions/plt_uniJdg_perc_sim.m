function plt_uniJdg_perc_sim(spatialD, true_pC1, fits_pC1, sim_unityJdg,...
    SD_simResp_pC1, numT_AV, ttl, bool_save)
if nargin < 8; bool_save = 0; end
if nargin < 7; ttl = '';end

global lenC lenP lenD
x_bds          = [spatialD(1) - 2, spatialD(end) + 2]; 
y_bds          = [-0.05,1.05];
x_ticks_AVdiff = spatialD(1:2:end); 
y_ticks        = 0:0.25:1;
lw             = 5; %lineWidth
fs_lbls        = 25;
fs_lgds        = 20;
jitter_prepost = [0.5, -0.5];
cMap_unity     = [0.65, 0.65,0.65;0.1,0.1,0.1];
lgd_cond_phase = {'Pre-learning', 'Post-learning'};
lgd_pos        = [0.35 0.15 0.05 0.1; 0.80 0.15 0.05 0.1]; 

figure
for i = 1:lenC
    subplot(1,2,i)
    addBackground(x_bds, y_bds, x_ticks_AVdiff, y_ticks)
    for j = 1:lenP
        if ~isempty(sim_unityJdg)
            for m = 1:lenD
                plt(j) = errorbar(spatialD(m)+jitter_prepost(j), ...
                    sim_unityJdg(i,j,m), SD_simResp_pC1(i,j,m),'-o',...
                    'MarkerSize',sqrt(numT_AV(m)).*1.5,'Color',cMap_unity(j,:),...
                    'MarkerFaceColor',cMap_unity(j,:),'MarkerEdgeColor',...
                    cMap_unity(j,:),'lineWidth',lw); hold on
            end
        end
        if ~isempty(true_pC1)
            plt(j) = plot(spatialD+jitter_prepost(j), squeeze(true_pC1(i,j,:)),...
                '-','lineWidth',lw,'Color',cMap_unity(j,:)); hold on;
        end
        if ~isempty(fits_pC1)
            plot(spatialD+jitter_prepost(j), squeeze(fits_pC1(i,j,:)),...
                '--','lineWidth',lw,'Color',cMap_unity(j,:)); hold on;
        end
    end
    %add legends
    legend([plt(1) plt(2)], lgd_cond_phase,'Position', lgd_pos(i,:),...
        'FontSize',fs_lgds); legend boxoff;
    xticks(x_ticks_AVdiff); xlim(x_bds); xlabel('Spatial discrepancy (V - A, dva)'); 
    ylabel(sprintf(['The probability of \nreporting ''common cause'''])); 
    yticks(y_ticks);ylim(y_bds);
    set(gca,'FontSize',fs_lbls);
end
sgtitle(ttl, 'FontSize',fs_lgds);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.62, 0.5]);
set(gcf,'PaperUnits','centimeters','PaperSize',[50 25]);
if bool_save; saveas(gcf, ['UnityJdg_sim_',ttl], 'pdf'); end

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


