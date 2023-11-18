function plot_bestfitPC1_diffSigma(pC1_diff, sigma_diff, delta_AIC_overall,...
    bestM_sub, bool_save)
global lenC lenM nS subjIs
% %if plot bootstrapped CI, load data
% cMAP  = [125,203,151;89,191,182;78,135,163;79,79,119;93,46,88;...
%          148,48,81;195,56,79;243,115,83;245,159,77;249,204,84;...
%          237,225,121;187,216,123;210,105,30;100,100,100;200,200,200]./255;
% cMAP_uni = {[242,182,126;100,195,187]./255,[215,133,76;91,144,169]./255};

%--------------------------------------------------------------------------
%sigma
x_bds_abs = ceil(max(abs(sigma_diff(:))));
x_bds   = [-x_bds_abs-1, x_bds_abs + 1];
x_ticks = [-x_bds_abs, 0, x_bds_abs]; 
x_lbls  = {{'\sigma_{AV,A,cong} - \sigma_{AV,A,prepost}',...
            '\sigma_{AV,A,incong} - \sigma_{AV,A,prepost}'}...
           {'\sigma_{AV,V,cong} - \sigma_{AV,V,prepost}',...
            '\sigma_{AV,V,incong} - \sigma_{AV,V,prepost}'}};
y_bds   = [min(pC1_diff(:))-0.2, max(pC1_diff(:)) + 0.2];
y_ticks = -0.25:0.25:0.75;
figure
for s = 1:lenM
    for i = 1:lenC
        subplot(lenM,lenC,(s-1)*lenM+i)
        addBackground(x_bds, y_bds, [x_bds(1), x_ticks, x_bds(end)],...
            [y_bds(1), y_ticks, y_bds(end)]);
        plot(x_bds, [0, 0], '-', 'lineWidth', 2,'Color',ones(1,3).*0.5); hold on;
        plot([0,0], y_bds, '-', 'lineWidth', 2,'Color',ones(1,3).*0.5); hold on;
        %select participants whose delta_AIC is greater than 7
        %(i.e., the complex model is much better than the simple model)
        selected_idx = [];
        for j = 1:nS
            delta_AIC = squeeze([min(delta_AIC_overall(j,1,:)), min(delta_AIC_overall(j,2,:))]);
            if diff(delta_AIC) < -7
                selected_idx = [selected_idx, j];
                h2 = scatter(sigma_diff(j,s,i), pC1_diff(j,i), 100, 'k','o',...
                    'MarkerFaceColor','r','MarkerFaceAlpha',0.5,'lineWidth',1); hold on
                h3 = text(sigma_diff(j,s,i)+0.2,pC1_diff(j,i),subjIs{j},'fontSize',10,...
                    'Color','r'); hold on
            else 
                h1 = scatter(sigma_diff(j,s,i), pC1_diff(j,i), 100, 'k','o',...
                    'MarkerFaceColor','k','MarkerFaceAlpha',0.5,'lineWidth',1); hold on
            end
        end
        
        %calculate mean delta sigma
        mean_sigma_diff = round(mean(sigma_diff(:,s,i)),2);
        SEM_sigma_diff  = round(std(sigma_diff(:,s,i))./sqrt(nS),2);
        
        mean_sigma_diff_selected = round(mean(sigma_diff(selected_idx,s,i)),2);
        SEM_sigma_diff_selected  = round(std(sigma_diff(selected_idx,s,i))./...
            sqrt(length(selected_idx)),2);
        
        xlim(x_bds); ylim(y_bds); yticks(y_ticks); xticks(x_ticks);
        xlabel(x_lbls{s}{i});
        if i == 1; ylabel('p_{C=1,post}-p_{C=1,pre}'); end
        legend([h1, h2], {['mean \Delta_{\sigma} = ', num2str(mean_sigma_diff),...
            ', SEM = ', num2str(SEM_sigma_diff)], ['mean \Delta_{\sigma} = ', ...
            num2str(mean_sigma_diff_selected), ', SEM = ', num2str(SEM_sigma_diff_selected)]},...
            'Location','northwest','fontSize',18); legend boxoff
        set(gca,'FontSize',25);
    end
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.65, 0.8]);
set(gcf,'PaperUnits','centimeters','PaperSize',[45 35]);
if bool_save; saveas(gcf, 'deltaSigma', 'pdf'); end


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
