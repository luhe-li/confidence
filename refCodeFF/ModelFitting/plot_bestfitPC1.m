function plot_bestfitPC1(pC1_bestModel, pC1_CI_lb, pC1_CI_ub,pC1_diff,...
    pC1_diff_CI_lb, pC1_diff_CI_ub, choice_c, subjIs, bool_save, idx_slc)
if nargin < 10; idx_slc = []; end
%if plot bootstrapped CI, load data
if isempty(pC1_CI_lb) || isempty(pC1_CI_lb); plt_CI_btst = 0;
else; plt_CI_btst = 1; end
x_lbl = {'pre-association', 'pre-dissociation'};
y_lbl = {'post-association', 'post-dissociation'};
cMAP  = [125,203,151;  89,191,182;  78,135,163;   79,79,119;    93,46,88;...
           148,48,81;   195,56,79;  243,115,83;  245,159,77;  249,204,84;...
         237,225,121; 187,216,123;  210,105,30; 100,100,100; 200,200,200;...
         110,200,250;  138,82,192]./255;
cMAP_uni = {[242,182,126;100,195,187]./255,[215,133,76;91,144,169]./255};
cond  = {'Congruent','Incongruent'};
x_bds = [0,1]; y_bds = [0,1]; x_ticks = 0:0.25:1; y_ticks = x_ticks;

%plot the two conditions separately
figure
for i = 1:length(cond)
    subplot(length(cond),1,i)
    addBackground(x_bds, y_bds, x_ticks, y_ticks)
    plot([0,1],[0,1],'k--','lineWidth',2); hold on
    for j = 1:size(pC1_bestModel,1)
        if plt_CI_btst
            errorbar(pC1_bestModel(j,i,1), pC1_bestModel(j,i,2),...
                pC1_bestModel(j,i,2) - pC1_CI_lb(j,i,2), ...
                pC1_CI_ub(j,i,2)-pC1_bestModel(j,i,2),...
                pC1_bestModel(j,i,1) - pC1_CI_lb(j,i,1), ...
                pC1_CI_ub(j,i,1)-pC1_bestModel(j,i,1),'-',...
                'MarkerFaceColor',cMAP(j,:),'MarkerEdgeColor',cMAP(j,:),...
                'lineWidth',2,'Color',cMAP(j,:)); hold on
        end
        scatter(pC1_bestModel(j,i,1), pC1_bestModel(j,i,2), 500, cMAP(j,:),...
            'd','MarkerFaceColor',cMAP(j,:),'MarkerFaceAlpha', 0.7,...
            'lineWidth',2); hold on
    end
    hold off; box off; axis equal; axis square;
    xticks(x_ticks); xlabel(sprintf(['Best-fitting p_{C=1} (',x_lbl{i},')']));
    yticks(y_ticks); ylabel(sprintf(['Best-fitting p_{C=1} (',y_lbl{i},')']));
    set(gca,'FontSize',18);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 1]);
set(gcf,'PaperUnits','centimeters','PaperSize',[20 32]);
if bool_save; saveas(gcf, 'pC1_bestModel_fig1', 'pdf'); end

%plot the two conditions together
jitter  = linspace(-0.25,0.25,size(pC1_diff,1));
x_bds   = [-1.5, 1.5];
x_ticks = -1:1;
y_bds   = [-0.45, 0.7];
y_ticks = -0.3:0.3:y_bds(end);
figure
addBackground(x_bds, y_bds, [x_bds(1), x_ticks, x_bds(end)], [y_bds(1), y_ticks, y_bds(end)]);
plot(x_bds, [0,0],'k--','lineWidth',2); hold on;
plot([0,0],y_bds,'k-','lineWidth',1.5); hold on;
for i = 1:size(pC1_diff,1)
    if strcmp(choice_c, 'rainbow'); cMAP_i = cMAP(idx_slc(i),:); 
    else
        if pC1_diff_CI_lb(i,1) > 0; cMAP_i = cMAP_uni{1}(1,:); 
        elseif pC1_diff_CI_ub(i,1) < 0; cMAP_i = cMAP_uni{1}(2,:);
        else; cMAP_i = ones(1,3).*0.5;
        end
    end
    errorbar(-1 + jitter(i), pC1_diff(i,1), pC1_diff(i,1) - pC1_diff_CI_lb(i,1),...
        pC1_diff_CI_ub(i,1) - pC1_diff(i,1),0,0,'-o','MarkerFaceColor','k',...
        'MarkerEdgeColor','k','lineWidth',2,'Color','k'); hold on
    scatter(-1 + jitter(i), pC1_diff(i,1), 500, 'k','d','MarkerFaceColor',...
        cMAP_i,'MarkerFaceAlpha', 0.6,'lineWidth',1); hold on
%     if strcmp(choice_c, 'rainbow')
%         plot([0,-1 + jitter(i)],[0,pC1_diff(i,1)],'Color', [cMAP_i, 0.5],...
%             'lineWidth',2.5); hold on;  
%     else
%         plot([0,-1 + jitter(i)],[0,pC1_diff(i,1)],'Color', [0,0,0,0.3],...
%             'lineWidth',2.5); hold on;  
%     end
%     text(-1 + jitter(i) - 0.03, pC1_diff(i,1), subjIs{i},'fontSize',10); hold on
    
    if strcmp(choice_c, 'rainbow'); cMAP_i = cMAP(idx_slc(i),:); 
    else
        if pC1_diff_CI_lb(i,2) > 0; cMAP_i = cMAP_uni{2}(1,:); 
        elseif pC1_diff_CI_ub(i,2) < 0; cMAP_i = cMAP_uni{2}(2,:);
        else; cMAP_i = ones(1,3).*0.5;
        end
    end
    errorbar(1 + jitter(i), pC1_diff(i,2), pC1_diff(i,2) - pC1_diff_CI_lb(i,2),...
        pC1_diff_CI_ub(i,2) - pC1_diff(i,2),0,0,'-o','MarkerFaceColor','k',...
        'MarkerEdgeColor','k','lineWidth',2,'Color','k'); hold on    
    scatter(1 + jitter(i), pC1_diff(i,2), 500, 'k','d','MarkerFaceColor',...
        cMAP_i,'MarkerFaceAlpha', 0.6,'lineWidth',1); hold on
%     if strcmp(choice_c, 'rainbow')
%         plot([0,1 + jitter(i)],[0,pC1_diff(i,2)],'Color', [cMAP_i, 0.5],...
%             'lineWidth',2.5); hold on;  
%     else
%         plot([0,1 + jitter(i)],[0,pC1_diff(i,2)],'Color', [0,0,0,0.3],...
%             'lineWidth',2.5); hold on;  
%     end
    if strcmp(choice_c, 'rainbow')
        plot([-1,1] + jitter(i), pC1_diff(i,:), 'Color', [cMAP_i, 0.5],...
            'lineWidth',2.5); hold on;  
    else
        plot([-1,1] + jitter(i), pC1_diff(i,:), 'Color', [0,0,0,0.3],...
            'lineWidth',2.5); hold on;          
    end
%     text(1 + jitter(i) - 0.03, pC1_diff(i,2), subjIs{i}, 'fontSize',10); hold on
end
xlim(x_bds); ylim(y_bds); xticks(x_ticks); yticks(y_ticks); 
ylabel(sprintf('Change of the common-cause prior\nafter a learning phase'));set(gca,'FontSize',25);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.65]);
set(gcf,'PaperUnits','centimeters','PaperSize',[35 35]);
if bool_save; saveas(gcf, ['pC1_bestModel_fig2_groups', subjIs{1}], 'pdf'); end


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
