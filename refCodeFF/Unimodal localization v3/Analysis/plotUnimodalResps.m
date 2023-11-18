% This script plots histograms for the unimodal localization responses
clear all; close all; clc; rng(1);

%% load and organize data
subjN    = 25; %[3,4,5,6,8,9,11,12,13,18,19,20,21,22,23,24,25]; %outliers: 15,16,17
subjI    = 'HHL'; %{'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','ZL','RE',...
%'MM','LL','CS','JH','MD','HHL'};%outliers: 'AD','SM','SX'
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'Unimodal localization v3/Data/',subjI]));
C        = load(['Unimodal_localization_sub', num2str(subjN), '.mat'],...
            'Unimodal_localization_data');
Adata     = C.Unimodal_localization_data{end}.data;
Vdata     = C.Unimodal_localization_data{end-1}.data;
Aloc_deg  = C.Unimodal_localization_data{end}.Distance; 
Vloc_deg  = C.Unimodal_localization_data{end-1}.initialDistance; 
numSloc   = length(Aloc_deg);
numTrials = size(Adata,2)/length(numSloc); %same for A and V = 30

%initialize matrices
numBtst = 1e3;
%mean localization responses for each stimulus location
[meanLocResp_A, meanLocResp_V] = deal(NaN(1, numSloc));
%mean localization responses for each bootstrapped dataset
[meanLocResp_A_btst, meanLocResp_V_btst] = deal(NaN(numBtst, numSloc));
%confidence intervals
[meanLocResp_A_95CI, meanLocResp_V_95CI] = deal(NaN(2, numSloc));
idx_95CI_lb = floor(0.025*numBtst);
idx_95CI_ub = ceil(0.975*numBtst);

for i = 1:numSloc
    idx_A = find(abs(Adata(1,:) - Aloc_deg(i)) < 1e-3);
    idx_V = find(abs(Vdata(1,:) - Vloc_deg(i)) < 1e-3);
    meanLocResp_A(i) = mean(Adata(2, idx_A));
    meanLocResp_V(i) = mean(Vdata(3, idx_V));
    for j = 1:numBtst
        idx_A_btst  = idx_A(randi([1 numTrials/numSloc],[1 numTrials/numSloc]));
        idx_V_btst  = idx_V(randi([1 numTrials/numSloc],[1 numTrials/numSloc]));
        AlocResp_btst = Adata(2, idx_A_btst);
        VlocResp_btst = Vdata(3, idx_V_btst);
        meanLocResp_A_btst(j,i) = mean(AlocResp_btst);
        meanLocResp_V_btst(j,i) = mean(VlocResp_btst);
    end
    meanLocResp_A_sort = sort(meanLocResp_A_btst(:,i));
    meanLocResp_V_sort = sort(meanLocResp_V_btst(:,i));
    meanLocResp_A_95CI(:,i) = meanLocResp_A_sort([idx_95CI_lb, idx_95CI_ub]);
    meanLocResp_V_95CI(:,i) = meanLocResp_V_sort([idx_95CI_lb, idx_95CI_ub]);
end

%% fit a linear regression
reg_A      = polyfit(Vloc_deg, meanLocResp_A, 1);
reg_V      = polyfit(Vloc_deg, meanLocResp_V, 1);
%for each bootstrapped dataset, compute a linear regression
reg_A_btst = arrayfun(@(idx) polyfit(Vloc_deg, meanLocResp_A_btst(idx,:), 1)',...
                1:numBtst, 'UniformOutput', false);
reg_V_btst = arrayfun(@(idx) polyfit(Vloc_deg, meanLocResp_V_btst(idx,:), 1)',...
                1:numBtst, 'UniformOutput', false);
%convert cells to matrices
reg_A_btst = cell2mat(reg_A_btst);
reg_V_btst = cell2mat(reg_V_btst);

%find 95% confidence intervals for the slope and the intercept
reg_A_slope_sort = sort(reg_A_btst(1,:)); 
reg_A_slope_95CI = reg_A_slope_sort([idx_95CI_lb, idx_95CI_ub]);
reg_A_intc_sort  = sort(reg_A_btst(2,:)); 
reg_A_intc_95CI  = reg_A_intc_sort([idx_95CI_lb, idx_95CI_ub]);

reg_V_slope_sort = sort(reg_V_btst(1,:)); 
reg_V_slope_95CI = reg_V_slope_sort([idx_95CI_lb, idx_95CI_ub]);
reg_V_intc_sort  = sort(reg_V_btst(2,:)); 
reg_V_intc_95CI  = reg_V_intc_sort([idx_95CI_lb, idx_95CI_ub]);

%% plot the figure
cMAP     = [255,178,178;255,127,127;229,0,0;153,0,0]./255;
x_max    = max(abs(Adata(2,:))); x_bds = [-x_max-5, x_max + 5];
y_bds    = [0,25]; x_ticks = [x_bds(1), Vloc_deg, x_bds(end)];
y_ticks  = 0:10:20;
modality = {'A','V'};
Modality = {'Auditory', 'Visual'};

figure
for m = 1:2 %modalities
    subplot(1,2,m)
    addBackground(x_bds, y_bds, x_ticks, y_ticks)
    for i = 1:numSloc
        eval(['d = ', modality{m},'data;']);
        eval(['loc = ', modality{m},'loc_deg(i);']);
        eval(['mR = meanLocResp_', modality{m},'(i);']);
        idx = abs(d(1,:) - loc) < 1e-3;
        histogram(d(end - 1,idx), -x_max:1:x_max,'FaceColor',cMAP(i,:),...
            'FaceAlpha', 0.3, 'EdgeColor',cMAP(i,:)); hold on
        if i == 1
            plot([Vloc_deg(i), Vloc_deg(i)], y_bds, 'Color',cMAP(i,:),...
                'lineWidth', 2, 'lineStyle', '-'); 
        else
            plot([Vloc_deg(end), Vloc_deg(end)], y_bds, 'Color', cMAP(i,:),...
                'lineWidth', 2, 'lineStyle', '-'); 
        end
        plot([mR,mR], [0, y_bds(end)], 'Color',...
            cMAP(i,:), 'lineWidth', 2, 'lineStyle', '--'); 
        text(x_bds(1), y_bds(end)-1, subjI, 'fontSize', 18); 
    end
    hold off; box off; xticks(Vloc_deg); xlim(x_bds); yticks(y_ticks);
    xlabel([Modality{m},' localization (dva)']);
    ylim(y_bds); ylabel('Counts'); set(gca,'FontSize',20);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.5]);
set(gcf,'PaperUnits','centimeters','PaperSize',[35 20]);
%saveas(gcf, ['UnimodalLocResps_S', num2str(subjN)], 'pdf'); 

%% plot the fitted linear regression
%confidence interval
polyval_btst_A = arrayfun(@(idx) polyval(reg_A_btst(:,idx)', Vloc_deg, 1)',...
    1:numBtst, 'UniformOutput', false);
polyval_btst_A = cell2mat(polyval_btst_A);
max_polyval_btst_A = max(polyval_btst_A');
min_polyval_btst_A = min(polyval_btst_A');

polyval_btst_V = arrayfun(@(idx) polyval(reg_V_btst(:,idx)', Vloc_deg, 1)',...
    1:numBtst, 'UniformOutput', false);
polyval_btst_V = cell2mat(polyval_btst_V);
max_polyval_btst_V = max(polyval_btst_V');
min_polyval_btst_V = min(polyval_btst_V');

axes_lb = min([round(meanLocResp_A(1)),Vloc_deg(1)]) - 4;
axes_ub = max([round(meanLocResp_A(end)),Vloc_deg(end)]) + 4;

figure
for m = 1:2
    subplot(1,2,m)
    addBackground([axes_lb, axes_ub], [axes_lb, axes_ub],...
        [axes_lb, Vloc_deg,axes_ub],[axes_lb, Vloc_deg, axes_ub]);
    %identity line
    plot(Vloc_deg, Vloc_deg, 'k--'); hold on
    %confidence interval
    eval(['min_polyval_btst = min_polyval_btst_', modality{m}, ';']);
    eval(['max_polyval_btst = max_polyval_btst_', modality{m}, ';']);
    patch([Vloc_deg, fliplr(Vloc_deg)], [min_polyval_btst, ...
        fliplr(max_polyval_btst)], [0.5, 0.5,0.5],'FaceAlpha',0.3,...
        'EdgeColor','none'); 
    plot(Vloc_deg, polyval(eval(['reg_', modality{m}]), Vloc_deg),...
        'Color','k','lineWidth',2); 
    for i = 1:numSloc
        eval(['mR = meanLocResp_', modality{m},'(i);']);
        eval(['mR_CI_lb = meanLocResp_', modality{m},'_95CI(1,i);']);
        eval(['mR_CI_ub = meanLocResp_', modality{m},'_95CI(2,i);']);
        errorbar(Vloc_deg(i), mR, mR -mR_CI_lb, mR_CI_ub - mR,'s',...
            'MarkerFaceColor', cMAP(i,:), 'MarkerEdgeColor', cMAP(i,:),...
            'Color',cMAP(i,:),'lineWidth',3,'MarkerSize',10);
    end
    if m == 1;text(axes_lb+.5, axes_ub-1, subjI, 'fontSize', 15); end
    hold off; axis equal; grid on; box off; xticks(Vloc_deg); 
    xlim([axes_lb, axes_ub]); ylim([axes_lb, axes_ub]); yticks(Vloc_deg);
    xlabel('Visual stimulus location (dva)');
    ylabel(['Mean ', Modality{m},' localization (dva)']); 
    set(gca,'FontSize',20);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.5]);
set(gcf,'PaperUnits','centimeters','PaperSize',[35 20]);
saveas(gcf, ['UnimodalLocResps_linearReg_S', num2str(subjN)], 'pdf'); 

%% save the data
data_summary = {[subjN, subjI], {meanLocResp_A, meanLocResp_V},...
    {meanLocResp_A_btst, meanLocResp_V_btst}, {meanLocResp_A_95CI, ...
    meanLocResp_V_95CI}, {reg_A, reg_V}, {reg_A_btst, reg_V_btst}, ...
    {reg_A_slope_95CI, reg_V_slope_95CI}, {reg_A_intc_95CI, reg_V_intc_95CI}};
save(['UnimodalLocalization_dataSummary_S', num2str(subjN)],'data_summary');



