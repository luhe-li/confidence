clear all; close all; clc
subjNs     = [   3,   4,   5,   6,   8,   9,  11,  12,  13,  18,  19,  20,...
                21,  22,  23,  24,  25]; %15AD, 16SM, 17SX are identified as outliers
subjIs     = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','ZL','RE','MM',...
              'LL','CS','JH','MD','HHL'};
nS         = length(subjNs);
spatialD   = -24:8:24;
numT_AV    = [20,40,60,80,60,40,20];
lenD       = length(spatialD);
cond       = {'congruent','incongruent'};
lenC       = length(cond);
phase      = {'pre','post'}; 
lenP       = length(phase);
modality   = {'A','V'};
lenM       = length(modality);

unityJdg      = NaN(nS, lenC, lenP, lenD);
diff_unityJdg = NaN(nS, lenC, lenD);
VE            = NaN(nS, lenC, lenP, lenM, lenD);
diff_VE       = NaN(nS, lenC, lenM, lenD);

for i = 1:nS
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'Bimodal localization v2/Data/', subjIs{i}]));
    C        = load(['Summary_VE_sub', num2str(subjNs(i)), '.mat'], 'Summary_VE');
    %size: 2 conditions x 2 phases x 7 spatial discrepancies
    unityJdg(i,:,:,:)    = C.Summary_VE{8}; 
    diff_unityJdg(i,:,:) = squeeze(unityJdg(i,:,2,:) - unityJdg(i,:,1,:));
    %size: 2 conditions x 2 phases x 2 modalities x 7 spatial discrepancies
    VE(i,:,:,:,:)        = C.Summary_VE{2}; 
    diff_VE(i,:,:,:)     = squeeze(VE(i,:,2,:,:) - VE(i,:,1,:,:)); %only auditory localizations are selected
end

%% plot the grouped VE
VE_Ashift      = squeeze(VE(:,:,:,1,:));
VE_Vshift      = squeeze(VE(:,:,:,2,:));
VE_Ashift_grp  = squeeze(mean(VE_Ashift, 1));
VE_Ashift_SEM  = squeeze(std(VE_Ashift,1))/sqrt(nS);
VE_Vshift_grp  = squeeze(mean(VE_Vshift, 1));
VE_Vshift_SEM  = squeeze(std(VE_Ashift,1))/sqrt(nS);
unityJdg_grp   = squeeze(mean(unityJdg,1));
unityJdg_SEM   = unityJdg_grp.*(1-unityJdg_grp)./sqrt(nS);

%combine across conditions and phase
VE_Ashift_squeeze_temp = squeeze(mean(mean(VE_Ashift,2),3));
VE_Ashift_squeeze = VE_Ashift_squeeze_temp(:, 4:end);
VE_Ashift_squeeze(:,2:end) = VE_Ashift_squeeze(:,2:end)./2 + (-fliplr(VE_Ashift_squeeze_temp(:,1:3))./2); 
VE_Ashift_squeeze_grp = squeeze(mean(VE_Ashift_squeeze, 1));
VE_Ashift_squeeze_SEM = squeeze(std(VE_Ashift_squeeze,1))/sqrt(nS);

VE_Vshift_squeeze_temp = squeeze(mean(mean(VE_Vshift,2),3));
VE_Vshift_squeeze = VE_Vshift_squeeze_temp(:, 4:end);
VE_Vshift_squeeze(:,2:end) = -VE_Vshift_squeeze(:,2:end)./2 + fliplr(VE_Vshift_squeeze_temp(:,1:3))./2; 
VE_Vshift_squeeze_grp = squeeze(mean(VE_Vshift_squeeze, 1));
VE_Vshift_squeeze_SEM = squeeze(std(VE_Vshift_squeeze,1))/sqrt(nS);

unityJdg_squeeze_temp = squeeze(mean(mean(unityJdg,2),3));
unityJdg_squeeze = unityJdg_squeeze_temp(:, 4:end);
unityJdg_squeeze(:,2:end) = unityJdg_squeeze(:,2:end)./2 + fliplr(unityJdg_squeeze_temp(:,1:3))./2; 
unityJdg_squeeze_grp = squeeze(mean(unityJdg_squeeze, 1));
unityJdg_squeeze_SEM = unityJdg_squeeze_grp.*(1-unityJdg_squeeze_grp)./sqrt(nS);

%% Plot grouped VE & proportion of reporting C = 1 
y_bds          = [-0.05,1.05];
x_ticks_AVdiff = spatialD; 
y_ticks        = 0:0.25:1;
y_bds_VE       = [spatialD(1)-2, spatialD(end)+2];
y_ticks_VE     = [spatialD(find(y_bds_VE(1) < spatialD,1):find(y_bds_VE(end) > spatialD,1,'last'))];
x_bds_VE       = [spatialD(1)-3, spatialD(end)+3];
x_ticks_VE     = [x_bds_VE(1), spatialD, x_bds_VE(end)];
cb             = [65,105,225]./255; cr = [0.85, 0.23, 0.28];
cMap_VE        = {min(cb.*2,1),cb; [247,191,190]./255, cr};
cMap_unity     = [0.65, 0.65,0.65;0.1,0.1,0.1];
lgd_pos        = [0.3 0.8 0.05 0.1; 0.74 0.8 0.05 0.1];
jitter_prepost = [-2, -0.67;0.67, 2];
fs             = 16.5;
marker_cond    = {'diamond','square'};
y_lbl          = {'Auditory ventriloquism effects (deg)',...
                  'Visual ventriloquism effects (deg)'};

figure
subplot(1,3,1) %unity judgment
addBackground(x_bds_VE, y_bds, x_ticks_VE, [y_bds(1),y_ticks,y_bds(end)])
for i = 1:lenC
    for j = 1:lenP
        for k = 1:length(spatialD)
            h(i,j) = errorbar(spatialD(k)+jitter_prepost(i,j), ...
                unityJdg_grp(i,j,k), unityJdg_SEM(i,j,k),'Marker',marker_cond{i},...
                'MarkerSize',sqrt(numT_AV(k)).*2,'Color',cMap_unity(j,:),...
                'MarkerFaceColor',cMap_unity(j,:),'MarkerEdgeColor',...
                cMap_unity(j,:),'lineWidth',2); hold on
        end
        plot(spatialD+jitter_prepost(i,j), squeeze(unityJdg_grp(i,j,:)),...
                '-','lineWidth',2,'Color',cMap_unity(j,:)); hold on;
    end
end
%add legends
text(x_bds_VE(1) + 0.2, 1.02, ['Group data (N=',num2str(nS),')'],...
    'FontSize', fs); hold on;
% legend([h(1,1), h(1,2), h(2,1),h(2,2)], {'Pre-association', 'Post-association',...
%     'Pre-dissociation', 'Post-dissociation'},'FontSize',20); legend boxoff;
xticks(spatialD); xlim(x_bds_VE); xlabel('Spatial discrepancy (V - A, deg)'); 
ylabel(sprintf(['The probability \n of reporting a common cause'])); 
yticks(y_ticks);ylim(y_bds); set(gca,'FontSize',fs);

subplot(1,3,2) %A
addBackground(x_bds_VE, y_bds_VE, x_ticks_VE, [y_bds_VE(1), y_ticks_VE,y_bds_VE(end)])
plot([spatialD(1),spatialD(end)], [spatialD(1),spatialD(end)],'k--',...
                    'Color',ones(1,3).*0.5,'lineWidth',2); hold on;
plot(spatialD, zeros(1,lenD), 'k:','Color', ones(1,3).*0.5,'lineWidth',2); 
for i = 1:lenC
    for j = 1:lenP
        for k = 1:length(spatialD)
            h(i,j) = errorbar(spatialD(k) + jitter_prepost(i,j), ...
                squeeze(VE_Ashift_grp(i,j,k)), squeeze(VE_Ashift_SEM(i,j,k)),...
                'Marker',marker_cond{i}, 'MarkerSize',sqrt(numT_AV(k)).*2,...
                'Color',cMap_VE{1,j}, 'MarkerFaceColor',cMap_VE{1,j},...
                'MarkerEdgeColor',cMap_VE{1,j}, 'lineWidth',2); hold on;
        end
        plot(spatialD + jitter_prepost(i,j), squeeze(VE_Ashift_grp(i,j,:)),...
            'lineWidth',2, 'Color',cMap_VE{1,j});
    end
end
%add legends
text(x_bds_VE(1) + 0.2, y_bds_VE(end)-1, ['Group data (N=',num2str(nS),')'],...
    'FontSize', fs); 
xticks(spatialD); xlim(x_bds_VE); xlabel('Spatial discrepancy (V-A, deg)');
yticks(y_ticks_VE); ylim(y_bds_VE); ylabel(sprintf(y_lbl{1})); 
set(gca,'FontSize',fs);

subplot(1,3,3) %V
addBackground(x_bds_VE, y_bds_VE, x_ticks_VE, [y_bds_VE(1), y_ticks_VE,y_bds_VE(end)])
plot([spatialD(1), spatialD(end)], [spatialD(end), spatialD(1)],'k--',...
     'Color',ones(1,3).*0.5,'lineWidth',2); hold on;
plot(spatialD, zeros(1,lenD), 'k:','Color', ones(1,3).*0.5,'lineWidth',2); 
for i = 1:lenC
    for j = 1:lenP
        for k = 1:length(spatialD)
            h(i,j) = errorbar(spatialD(k) + jitter_prepost(i,j), ...
                squeeze(VE_Vshift_grp(i,j,k)), squeeze(VE_Vshift_SEM(i,j,k)),...
                'Marker',marker_cond{i}, 'MarkerSize',sqrt(numT_AV(k)).*2,...
                'Color',cMap_VE{2,j}, 'MarkerFaceColor',cMap_VE{2,j},...
                'MarkerEdgeColor',cMap_VE{2,j}, 'lineWidth',2); hold on;
        end
        plot(spatialD + jitter_prepost(i,j), squeeze(VE_Vshift_grp(i,j,:)),...
            'lineWidth',2, 'Color',cMap_VE{2,j});
    end
end
%add legends
text(x_bds_VE(1) + 0.2, y_bds_VE(end)-1, ['Group data (N=',num2str(nS),')'],...
    'FontSize', fs); 
legend([h(1,1), h(1,2), h(2,1),h(2,2)], {'Congruent pre-learning',...
    'Congruent post-learning', 'Incongruent pre-learning', ...
    'Incongruent post-learning'},'Location','southeast',...
    'FontSize',fs); legend boxoff;
xticks(spatialD); xlim(x_bds_VE); xlabel('Spatial discrepancy (V-A, deg)');
yticks(y_ticks_VE); ylim(y_bds_VE); ylabel(sprintf(y_lbl{2})); 
set(gca,'FontSize',fs);

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.52]);
set(gcf,'PaperUnits','centimeters','PaperSize',[60 25]);
saveas(gcf, 'VE_UnityJdg_groupMean', 'pdf');

%% combine the conditions and phases
spatialD_abs    = unique(abs(spatialD));
x_bds_VE        = [spatialD_abs(1)-3, spatialD_abs(end)+3];
x_ticks_VE      = [x_bds_VE(1), spatialD_abs, x_bds_VE(end)];
numT_AV_squeeze = numT_AV(4:end) + [0,fliplr(numT_AV(1:3))];
cMap_VE_squeeze = [87,137,224; 250,115,125]./255;

figure
subplot(1,3,1) %unity judgment
addBackground(x_bds_VE, y_bds, x_ticks_VE, [y_bds(1), y_ticks, y_bds(end)])
for k = 1:length(spatialD_abs)
    errorbar(spatialD_abs(k), unityJdg_squeeze_grp(k), unityJdg_squeeze_SEM(k),...
        '-o','MarkerSize',sqrt(numT_AV_squeeze(k)).*2,'Color',ones(1,3).*0.4,...
        'MarkerFaceColor',ones(1,3).*0.4,'MarkerEdgeColor',ones(1,3).*0.4,...
        'lineWidth',2); hold on
end
plot(spatialD_abs, unityJdg_squeeze_grp, '-','lineWidth',2,'Color',...
    ones(1,3).*0.4); hold on;
%add legends
text(x_bds_VE(1) + 0.2, 1.02, ['Group data (N=',num2str(nS),')'],...
    'FontSize', fs); hold on;
xticks(spatialD_abs); xlim(x_bds_VE); xlabel('Absolute spatial discrepancy (deg)'); 
ylabel(sprintf(['The probability \n of reporting a common cause'])); 
yticks(y_ticks);ylim(y_bds); set(gca,'FontSize',fs);

subplot(1,3,2) %A
addBackground(x_bds_VE, x_bds_VE, x_ticks_VE, x_ticks_VE)
plot([spatialD_abs(1),spatialD_abs(end)], [spatialD_abs(1),spatialD_abs(end)],'k--',...
	'Color',ones(1,3).*0.5,'lineWidth',2); hold on;
plot([spatialD_abs(1),spatialD_abs(end)], zeros(1,2), 'k:','Color',...
    ones(1,3).*0.5,'lineWidth',2); 
for k = 1:length(spatialD_abs)
	errorbar(spatialD_abs(k), VE_Ashift_squeeze_grp(k), VE_Ashift_squeeze_SEM(k),...
        '-o','MarkerSize',sqrt(numT_AV_squeeze(k)).*2,...
        'Color',cMap_VE_squeeze(1,:), 'MarkerFaceColor',cMap_VE_squeeze(1,:),...
        'MarkerEdgeColor',cMap_VE_squeeze(1,:), 'lineWidth',2); hold on;
end
plot(spatialD_abs, VE_Ashift_squeeze_grp, 'lineWidth',2, 'Color',cMap_VE_squeeze(1,:));
%add legends
text(x_bds_VE(1) + 0.2, y_bds_VE(end)-1, ['Group data (N=',num2str(nS),')'],...
    'FontSize', fs); 
xticks(spatialD_abs); xlim(x_bds_VE); xlabel('Absolute spatial discrepancy (deg)'); 
yticks(spatialD_abs(1:2)); ylim([-3, 11]); ylabel(sprintf(y_lbl{1})); 
set(gca,'FontSize',fs);

subplot(1,3,3) %V
addBackground(x_bds_VE, x_bds_VE, x_ticks_VE, [-0.3, 0,0.8,1.1])
plot([spatialD_abs(1),spatialD_abs(end)], [spatialD_abs(1),spatialD_abs(end)],'k--',...
	'Color',ones(1,3).*0.5,'lineWidth',2.5); hold on;
plot([spatialD_abs(1),spatialD_abs(end)], zeros(1,2), 'k:','Color',...
    ones(1,3).*0.5,'lineWidth',2); 
for k = 1:length(spatialD_abs)
	errorbar(spatialD_abs(k), VE_Vshift_squeeze_grp(k), VE_Vshift_squeeze_SEM(k),...
        '-o','MarkerSize',sqrt(numT_AV_squeeze(k)).*2,'Color',cMap_VE_squeeze(2,:),...
        'MarkerFaceColor',cMap_VE_squeeze(2,:), 'MarkerEdgeColor',cMap_VE_squeeze(2,:), ...
        'lineWidth',2); hold on;
end
plot(spatialD_abs, VE_Vshift_squeeze_grp, 'lineWidth',2, 'Color',cMap_VE_squeeze(2,:));
%add legends
text(x_bds_VE(1) + 0.2, y_bds_VE(end)-1, ['Group data (N=',num2str(nS),')'],...
    'FontSize', fs); 
xticks(spatialD_abs); xlim(x_bds_VE); xlabel('Absolute spatial discrepancy (deg)'); 
yticks([0,0.8]); ylim([-0.3, 1.1]);  ylabel(sprintf(y_lbl{2})); 
set(gca,'FontSize',fs);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.52]);
set(gcf,'PaperUnits','centimeters','PaperSize',[60 25]);
saveas(gcf, 'VE_UnityJdg_groupMean_squeeze', 'pdf');

% %%
% figure
% for i = 1:lenC
%     for j = 1:lenP
%         unityJdg_ij = squeeze(unityJdg(:,i,j,:));
%         unityJdg_ij_squeeze = unityJdg_ij(:,4:end);
%         unityJdg_ij_squeeze(:,2:end) = unityJdg_ij_squeeze(:,2:end)./2 + fliplr(unityJdg_ij(:,1:3)./2);
%         unityJdg_ij_mean = mean(unityJdg_ij_squeeze,1);
%         unityJdg_ij_SEM  = unityJdg_ij_mean.*(1-unityJdg_ij_mean)./sqrt(nS);
%         errorbar(spatialD_abs+jitter_prepost(i,j), unityJdg_ij_mean, unityJdg_ij_SEM); hold on
%     end
% end
% legend({'Cong, pre-', 'Cong, post-', 'Incong, pre-', 'Incong, post-'}); 
% xticks(spatialD_abs);
% set(gca,'FontSize',20);

%% plot the correlation
%first exclude people who show no effect
addpath(genpath('/Users/hff/Desktop/NYU/Project2/Experiment code/ModelFitting'));
C = load('Results_modelComparison_overall.mat', 'ModelComparison');
bestM = C.ModelComparison{9};

diff_VE_auditory = squeeze(diff_VE(:,:,1,:));
[rho,pval] = corr(diff_unityJdg(:),diff_VE_auditory(:));
disp([rho,pval]);
slope = polyfit(diff_unityJdg(:),diff_VE_auditory(:),1);
cMAP  = [125,203,151;  89,191,182;  78,135,163;   79,79,119;    93,46,88;...
           148,48,81;   195,56,79;  243,115,83;  245,159,77;  249,204,84;...
         237,225,121; 187,216,123;  210,105,30; 100,100,100; 200,200,200;...
         110,200,250;  138,82,192]./255;
marker  = {'o','o'};%{'o','s'};
x_bds   = [-1,1].*abs(max(diff_unityJdg(:))) + [-0.1,0.1];
y_bds   = [-15,15];
x_ticks = [x_bds(1), x_bds(1)/2, 0, x_bds(end)/2, x_bds(end)];
y_ticks = linspace(y_bds(1), y_bds(end),5);
ii = 0;
     
figure
for i = 1:nS
    seg_idx  = find(bestM{i} == '-'); 
    bestM_subj_i = bestM{i}((seg_idx(1)+1):(seg_idx(2)-1));
    if strcmp(bestM_subj_i, 'fullModel') || strcmp(bestM_subj_i, 'samePC1pre')
        ii = ii + 1;
        for j = 1:lenC
            for k = 1:lenD
                scatter(diff_unityJdg(i,j,k), diff_VE_auditory(i,j,k), 100, marker{j},'filled',...
                    'MarkerFaceColor', cMAP(ii,:), 'MarkerEdgeColor', cMAP(ii,:).*0.5,...
                    'MarkerFaceAlpha',0.5); hold on
            end
        end
        diff_unityJdg_subi = diff_unityJdg(i,:,:);
        diff_VE_auditory_subi = diff_VE_auditory(i,:,:);
    end
end
plot(x_bds, polyval(slope,x_bds), 'lineWidth',3,'Color','k'); hold on
text(x_bds(1)+0.05, y_bds(2)-2, ['r = ',num2str(round(rho,3)),', p = <0.001'],...
    'fontSize',15); hold off; box off
xlim(x_bds); ylim(y_bds); xticks(x_ticks); yticks(y_ticks);
xlabel(sprintf('Change in the proportion \nof reporting a common cause'));
ylabel(sprintf('Change in the \nauditory localization shifts'));
set(gca,'FontSize',20);




