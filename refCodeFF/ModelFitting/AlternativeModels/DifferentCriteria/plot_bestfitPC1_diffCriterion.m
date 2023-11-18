function plot_bestfitPC1_diffCriterion(pC1_diff, criterion_diff, delta_AIC_overall,...
    bestM_sub, choice_c, bool_save)
global lenC nS subjIs
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
y_bds   = [-0.75, 0.95];
y_ticks = -0.75:0.25:y_bds(end);

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
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.5]);
set(gcf,'PaperUnits','centimeters','PaperSize',[25 20]);
if bool_save; saveas(gcf, 'pC1_diffCriterion', 'pdf'); end

%--------------------------------------------------------------------------
%sigma
x_bds_abs = 1;%ceil(max(abs(criterion_diff(:)))); 
x_bds   = [-x_bds_abs, x_bds_abs];
x_ticks = linspace(-x_bds_abs,x_bds_abs,5); 
x_lbls  = {'\epsilon_{post-cong} - \epsilon_{pre-cong}',...
            '\epsilon_{post-incong} - \epsilon_{pre-incong}'};
y_bds   = [-0.75, 0.75];
y_ticks = -0.5:0.5:0.5;
figure
for i = 1:lenC
    subplot(1,lenC,i)
    addBackground(x_bds, y_bds, [x_bds(1), x_ticks, x_bds(end)],...
        [y_bds(1), y_ticks, y_bds(end)]);
    plot(x_bds, [0, 0], '-', 'lineWidth', 2,'Color',ones(1,3).*0.5); hold on;
    plot([0,0], y_bds, '-', 'lineWidth', 2,'Color',ones(1,3).*0.5); hold on;
    %select participants whose delta_AIC is greater than 7
    %(i.e., the complex model is much better than the simple model)
    for j = 1:nS
        if j ~= 5
            scatter(criterion_diff(j,i), pC1_diff(j,i), 100, 'k','o',...
                'MarkerFaceColor','k','MarkerFaceAlpha',0.5,'lineWidth',1); hold on
            text(criterion_diff(j,i)+0.02,pC1_diff(j,i),subjIs{j},'fontSize',10,...
                        'Color','r'); hold on
        end
    end
    hold off;
    xlim(x_bds); ylim(y_bds); yticks(y_ticks); xticks(x_ticks);
    xlabel(x_lbls{i});
    if i == 1; ylabel('p_{C=1,post}-p_{C=1,pre}'); end
    set(gca,'FontSize',25);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.65, 0.4]);
set(gcf,'PaperUnits','centimeters','PaperSize',[45 20]);
if bool_save; saveas(gcf, 'deltaCriterion', 'pdf'); end


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
