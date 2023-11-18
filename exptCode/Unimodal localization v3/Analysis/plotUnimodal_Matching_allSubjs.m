%This script aims to compare the mean localization responses collected from
%the unimodal spatial localization task with the perceptually matched
%auditory locations measured from the bimodal spatial discrimination task
clear all; close all; clc; rng(1)
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'Unimodal localization v3/Data']));
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/Matching/Data/']));

%% load data files
subjNs        = [   3,   4,   5,   6,   8,   9,  11,  12,  13,  15,  16,...
                   17,  18,  19,  20,  21,  22,  23,  24,  25];
subjIs        = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','AD','SM',...
                 'SX','ZL','RE','MM','LL','CS','JH','MD','HHL'};
lenS          = length(subjNs);
Vloc          = -12:8:12;
nBtst         = 1000;
[Reg_A_uni, Reg_A_slope_95CI_uni,  Reg_A_intc_95CI_uni,...
    Reg_V_uni, Reg_V_slope_95CI_uni,  Reg_V_intc_95CI_uni] = deal(NaN(lenS, 2)); 

for i = 1:lenS
    B = load(['UnimodalLocalization_dataSummary_S', num2str(subjNs(i)),...
        '.mat'], 'data_summary');
    Reg_A_uni(i,:) = B.data_summary{end-3}{1};
    Reg_A_slope_95CI_uni(i,:) = B.data_summary{end-1}{1};
    Reg_A_intc_95CI_uni(i,:) = B.data_summary{end}{1};

    Reg_V_uni(i,:) = B.data_summary{end-3}{2};
    Reg_V_slope_95CI_uni(i,:) = B.data_summary{end-1}{2};
    Reg_V_intc_95CI_uni(i,:) = B.data_summary{end}{2};
end

%find outliers
outlierC               = 3; subj_outlier = [];
Reg_slope_absDiff      = abs(Reg_A_uni(:,1) - Reg_V_uni(:,1));
Reg_slope_mean_absDiff = mean(Reg_slope_absDiff);
Reg_slope_std_absDiff  = std(Reg_slope_absDiff);
subj_outlier           = find(Reg_slope_absDiff > outlierC*Reg_slope_std_absDiff | Reg_slope_absDiff < -outlierC*Reg_slope_std_absDiff);

Reg_intc_absDiff      = abs(Reg_A_uni(:,2) - Reg_V_uni(:,2));
Reg_intc_mean_absDiff = mean(Reg_intc_absDiff);
Reg_intc_std_absDiff  = std(Reg_intc_absDiff);
subj_outlier          = [subj_outlier; find(Reg_intc_absDiff > outlierC*Reg_intc_std_absDiff | Reg_intc_absDiff < -outlierC*Reg_intc_std_absDiff)];
disp('Outliers: ');
for i = 1:length(subj_outlier)
    disp([num2str(subjNs(subj_outlier(i))), subjIs{subj_outlier(i)}]);
end

%display mean absolute difference
fprintf('Mean absolute difference (slope): %.4f; (intercept): %.4f\n',...
    Reg_slope_mean_absDiff, Reg_intc_mean_absDiff);
fprintf('STD absolute difference (slope): %.4f; (intercept): %.4f\n',...
    Reg_slope_std_absDiff, Reg_intc_std_absDiff);

%stats with outliers removed
subjNs_rm                 = 1:lenS;
subjNs_rm(subj_outlier)   = [];
Reg_slope_absDiff_rm      = abs(Reg_A_uni(subjNs_rm,1) - Reg_V_uni(subjNs_rm,1));
Reg_slope_mean_absDiff_rm = mean(Reg_slope_absDiff_rm);
Reg_slope_std_absDiff_rm  = std(Reg_slope_absDiff_rm);
Reg_intc_absDiff_rm       = abs(Reg_A_uni(subjNs_rm,2) - Reg_V_uni(subjNs_rm,2));
Reg_intc_mean_absDiff_rm  = mean(Reg_intc_absDiff_rm);
Reg_intc_std_absDiff_rm   = std(Reg_intc_absDiff_rm);

%display mean absolute difference after removing outliers
disp('After removing outliers');
fprintf('Mean absolute difference (slope): %.4f; (intercept): %.4f\n',...
    Reg_slope_mean_absDiff_rm, Reg_intc_mean_absDiff_rm);
fprintf('STD absolute difference (slope): %.4f; (intercept): %.4f\n',...
    Reg_slope_std_absDiff_rm, Reg_intc_std_absDiff_rm);

%% plot
x_bds   = [0,2];%[0.5, 1.5];
y_bds   = [0,2];
x_ticks = 0:0.5:2;%0.5:0.25:1.5;
y_ticks = 0:0.5:2;
cMAP    = [222,110,250; 171,83,193;123,57,139]./255;
%mstyle = {'o','s','d','>','<','^','v'};

figure
subplot(1,2,1)
addBackground(y_bds, y_bds, x_ticks, y_ticks)
plot(x_bds, y_bds, 'k--'); hold on
for i = 1:lenS
    if ismember(i, subj_outlier); ccMAP = cMAP(i==subj_outlier,:);
    else; ccMAP = [0.5, 0.5, 0.5]; end
    errorbar(Reg_V_uni(i,1), Reg_A_uni(i,1),...
        Reg_V_uni(i,1) - Reg_V_slope_95CI_uni(i,1), Reg_V_slope_95CI_uni(i,2) - ...
        Reg_V_uni(i,1), Reg_A_uni(i,1) - Reg_A_slope_95CI_uni(i,1),...
        Reg_A_slope_95CI_uni(i,2) - Reg_A_uni(i,1),'.', 'Color', ccMAP,...
        'lineWidth',2,'MarkerSize',10,'CapSize',8); hold on
end
xticks(x_ticks); xlim(x_bds);  yticks(y_ticks); ylim(y_bds); 
xlabel(sprintf(['Slope of the linear regression to \n visual localization responses (dva)']));
ylabel(sprintf(['Slope of the linear regression to \n auditory localization responses (dva)']));
axis square; box off; grid on; set(gca,'FontSize',20);

x_bds   = [-8,4];%[-4,4];
y_bds   = [-8,4];
x_ticks = -8:4:4;%-4:2:4;
y_ticks = -8:4:4;

subplot(1,2,2)
addBackground(x_bds, y_bds, x_ticks, y_ticks)
plot(x_bds, x_bds, 'k--'); hold on
for i = 1:lenS
    if ismember(i, subj_outlier);  ccMAP = cMAP(i==subj_outlier,:); 
    else; ccMAP = [0.5, 0.5, 0.5]; end
    errorbar(Reg_V_uni(i,2), Reg_A_uni(i,2),...
        Reg_V_uni(i,2) - Reg_V_intc_95CI_uni(i,1), Reg_V_intc_95CI_uni(i,2) - ...
        Reg_V_uni(i,2), Reg_A_uni(i,2) - Reg_A_intc_95CI_uni(i,1),...
        Reg_A_intc_95CI_uni(i,2) - Reg_A_uni(i,2),'.', 'Color', ccMAP,...
        'lineWidth',2,'MarkerSize', 10,'CapSize',8); hold on
end
xticks(x_ticks); xlim(x_bds);  yticks(y_ticks); ylim(y_bds); grid on
xlabel(sprintf(['Intercept of the linear regression to \n visual localization responses (dva)']));
ylabel(sprintf(['Intercept of the linear regression to \n auditory localization responses (dva)']));
axis square; box off; grid on; set(gca,'FontSize', 20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.5]);
set(gcf,'PaperUnits','centimeters','PaperSize',[35 20]);
saveas(gcf, 'Reg_unimodal_allSubjs', 'pdf');







