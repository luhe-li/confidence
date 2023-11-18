function plot_bestfitPC1_diffSigma(pC1_diff, sigma_diff, delta_AIC_overall,...
    bestM_sub, choice_c, bool_save)
global lenC lenM nS subjIs
%if plot bootstrapped CI, load data
cMAP  = [125,203,151;89,191,182;78,135,163;79,79,119;93,46,88;...
         148,48,81;195,56,79;243,115,83;245,159,77;249,204,84;...
         237,225,121;187,216,123;210,105,30;100,100,100;200,200,200]./255;
cMAP_uni = {[242,182,126;100,195,187]./255,[215,133,76;91,144,169]./255};

%plot the two conditions together
jitter  = linspace(-0.25,0.25,size(pC1_diff,1));
x_offset= [-1,1];
x_bds   = [-1.5, 1.5];
x_ticks = -1:1;
y_bds   = [-0.6, 0.95];
y_ticks = -0.5:0.25:y_bds(end);
figure
addBackground(x_bds, y_bds, [x_bds(1), x_ticks, x_bds(end)],...
    [y_bds(1), y_ticks, y_bds(end)]);
plot(x_bds, [0,0],'k--','lineWidth',2); hold on;
plot([0,0],y_bds,'k-','lineWidth',1.5); hold on;
for i = 1:nS
    for j = 1:lenC
        if strcmp(choice_c, 'rainbow'); cMAP_i = cMAP(i,:); 
        else
            if pC1_diff(i,j) > 0; cMAP_i = cMAP_uni{j}(1,:); 
            elseif pC1_diff(i,j) < 0; cMAP_i = cMAP_uni{j}(2,:);
            else; cMAP_i = ones(1,3).*0.5;
            end
        end
        scatter(x_offset(j) + jitter(i), pC1_diff(i,j), 500, 'k',...
            'd','MarkerFaceColor',cMAP_i,'MarkerFaceAlpha', 0.6,...
            'lineWidth',1); hold on
        plot([0,x_offset(j) + jitter(i)],[0,pC1_diff(i,j)],'Color',...
            [0,0,0,0.3],'lineWidth',2.5); hold on;
        text(x_offset(j) + jitter(i) - 0.03, pC1_diff(i,j), subjIs{i},...
            'fontSize',10); hold on
    end 
end
xlim(x_bds); ylim(y_bds); xticks([]); yticks(y_ticks); 
ylabel('p_{C=1,post}-p_{C=1,pre}');set(gca,'FontSize',25);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.6]);
set(gcf,'PaperUnits','centimeters','PaperSize',[30 30]);
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.75]);
% set(gcf,'PaperUnits','centimeters','PaperSize',[28 30]);
if bool_save; saveas(gcf, 'pC1_diffSigma', 'pdf'); end

%--------------------------------------------------------------------------
%sigma
x_bds_abs = 8;%ceil(max(abs(sigma_diff(:)))); 
x_bds   = [-x_bds_abs-1, x_bds_abs + 1];
x_ticks = [-x_bds_abs, 0, x_bds_abs]; 
x_lbls  = {{'\sigma_{AV,A,post-cong} - \sigma_{AV,A,pre}',...
            '\sigma_{AV,A,post-incong} - \sigma_{AV,A,pre}'}...
           {'\sigma_{AV,V,post-cong} - \sigma_{AV,V,pre}',...
            '\sigma_{AV,V,post-incong} - \sigma_{AV,V,pre}'}};
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
            delta_AIC = squeeze(delta_AIC_overall(j,:,bestM_sub(j,2)));
            if diff(delta_AIC) < -4
                selected_idx = [selected_idx, j];
                h2= scatter(sigma_diff(j,s,i), pC1_diff(j,i), 100, 'k','o',...
                    'MarkerFaceColor','r','MarkerFaceAlpha',0.5,'lineWidth',1); hold on
                h3 = text(sigma_diff(j,s,i)+0.2,pC1_diff(j,i),subjIs{j},'fontSize',10,...
                    'Color','r'); hold on
            else 
                h1=scatter(sigma_diff(j,s,i), pC1_diff(j,i), 100, 'k','o',...
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
