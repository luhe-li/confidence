%% set ground truth
clear all; close all; clc
%subject number (outlier participants: 15, 16, 17)
subjN_dict = [   3,   4,   5,   6,   8,   9,  11,  12,  13,...
                18,  19,  20,  21,  22,  23,  24, 25];
%subject initial (outlier participants: 'AD', 'SM', 'SX'
subjI_dict = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR',...
              'ZL','RE','MM','LL','CS','JH','MD','HHL'};
%the order of condition was counterbalanced across participants
%1: congruent; 2: incongruent
order_dict = [ 2,1; 1,2; 2,1; 1,2; 2,1; 1,2; 2,1; 2,1; 1,2;...
               2,1; 1,2; 2,1; 2,1; 1,2; 2,1; 1,2; 2,1];
           
%subject info
subjN             = 12; 
subjI             = 'ZY';
subj_idx          = find(subjN == subjN_dict);
sesNum            = order_dict(subj_idx,1); %congruent condition
          
%% load parameter estimate
%load model compariosn
addpath(genpath('/Users/hff/Desktop/NYU/Project2/Experiment code/ModelFitting/'));
B                 = load('Results_modelComparison_overall.mat', 'ModelComparison');
subjI_bestM       = B.ModelComparison{end};
%again load the best-fitting model
bestM             = subjI_bestM{subj_idx};
seg_idx           = find(bestM == '-'); 
model.M_pC1       = bestM((seg_idx(1)+1):((seg_idx(end)-1)));
if strcmp(model.M_pC1, 'fullModel') || strcmp(model.M_pC1, 'samePC1pre'); pC1_idx = [8,9];
else; pC1_idx = [8,8]; end
strategy          = bestM((seg_idx(end)+1):end);
seg_idx           = find(strategy == '_'); 
model.M_estimates = strategy((seg_idx(1)+1):((seg_idx(2)-1)));
model.M_unityJdg  = strategy((seg_idx(end)+1):end);

addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'ModelFitting/Fits/YZ/Fits_locResp_conditioned_unityJdg_jointly/']));
C = load(['ModelFitting_updatePrior_',model.M_pC1,'_', strategy,...
    '_sub', num2str(subjN), '.mat'], 'ModelFitting');
estimatedP        = C.ModelFitting{end}.P_f;
param.mu_P        = 0;
param.sigma_P     = 100;
param.a_A         = estimatedP(1);
param.b_A         = estimatedP(2);
param.pCommon     = estimatedP(pC1_idx(1));
param.c_unityJdg  = estimatedP(7);
alpha_pC1         = logspace(-2, -0.6, 5); 

addpath(genpath('/Users/hff/Desktop/NYU/Project2/Experiment code/Adaptation v2/Data'));
D = load('Analysis_learningPhase.mat', 'data');
data_learning    = D.data{3};
sigma_r          = D.data{4}(subj_idx);
sigma_A_learning = sqrt(data_learning.locVar{1}(subj_idx,1)^2 - sigma_r^2);
sigma_V_learning = sqrt(data_learning.locVar{2}(subj_idx,1)^2 - sigma_r^2);
if imag(sigma_A_learning)~=0; sigma_A_learning = 1e-2;end
if imag(sigma_V_learning)~=0; sigma_V_learning = 1e-2;end

sigma_A_learning_vec = linspace(.1, sigma_A_learning*0.5, 20);
sigma_V_learning_vec = linspace(1e-4, sigma_V_learning*1.2, 20);

%% stimulus pairs presented during the experiments
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Adaptation v2/Data/',subjI]));
F = load(['Adaptation_congruent_sub', num2str(subjN),'_session',...
    num2str(sesNum),'.mat'], 'Adaptation_data');
D.AVpairs(1,:)        = F.Adaptation_data{4}.arrangedLocs_deg; %A loc
D.AVpairs(2,:)        = F.Adaptation_data{3}.arrangedLocs_deg; %V loc
D.AVpairs(3,:)        = zeros(1, length(D.AVpairs(1,:)));
D.AVpairs(4,:)        = zeros(1, length(D.AVpairs(1,:)));
D.totalTrials         = size(D.AVpairs,2);
D.numSims             = 100;
D.bool_unity          = ~isnan(F.Adaptation_data{end}.unity);

%%
[propC1, pCommon_end] = deal(NaN(length(sigma_A_learning_vec), length(sigma_V_learning_vec),...
    length(alpha_pC1), length(D.numSims)));

for i = 1:length(sigma_A_learning_vec)
    disp(i)
    param.sigma_A_learning = sigma_A_learning_vec(i);
    for j= 1:length(sigma_V_learning_vec)
        param.sigma_V_learning = sigma_V_learning_vec(j);
        for k = 1:length(alpha_pC1)
            for l = 1:D.numSims
                [propC1(i,j,k,l), pCommon_end(i,j,k,l), ~] = simulate_learningPhase(...
                    param, alpha_pC1(k), D, model);
            end
        end
    end
end

%% plot
%load colormap
addpath(genpath('/Users/hff/Desktop/NYU/Project2/Experiment code/Simulations/cbrewer'));
mean_pCommon_end   = mean(pCommon_end,4);
simChange_pCommon  = mean_pCommon_end-param.pCommon;
empChange_pCommon  = estimatedP(pC1_idx(2)) - param.pCommon;
simReportingC1     = mean(propC1,3); 
empReportingC1     = sum(F.Adaptation_data{end}.unity==1)/...
                     sum(~isnan(F.Adaptation_data{end}.unity));

figure
min_v1   = min(mean_pCommon_end(:)); 
max_v1   = max(mean_pCommon_end(:));
cmap_cus = customize_colormap(max_v1, min_v1, param.pCommon); colormap(cmap_cus); 
RdBu     = flipud(cbrewer('div', 'RdBu', 51)); 
for i = 1:length(alpha_pC1)
    subplot(1,length(alpha_pC1),i)
    imagesc('XData',sigma_V_learning_vec,'YData',sigma_A_learning_vec,...
        'CData',squeeze(mean_pCommon_end(:,:,i))); 
    if estimatedP(11) > max_v1; caxis([min_v1, estimatedP(11)]);
    else; caxis([min_v1, max_v1]); end
    xlim([sigma_V_learning_vec(1), sigma_V_learning_vec(end)]);
    ylim([sigma_A_learning_vec(1), sigma_A_learning_vec(end)]);
    if i == 1
        yticks(round(sigma_A_learning_vec(1:6:end),1)); 
        ylabel('$\sigma_{AV,A}''$ (learning)','interpreter','latex'); 
    else; yticks([]);
    end
    xticks(round(sigma_V_learning_vec(1:6:end),2));
    title(['$\alpha_{p_{C=1}} = $', num2str(round(alpha_pC1(i),3))],...
        'interpreter','latex');
    if i == length(alpha_pC1); colorbar; caxis([min_v1, max_v1]); end
    xlabel('$\sigma_{AV,V}''$ (learning)','interpreter','latex');
    set(gca,'FontSize',15);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.40]);
set(gcf,'PaperUnits','centimeters','PaperSize',[70 25]);
saveas(gcf, ['Simtest_cong_', subjI], 'pdf'); 
    
%%
figure
max_v2 = max(simReportingC1(:)); min_v2 = min(simReportingC1(:));
cmap_cus2 = customize_colormap(max_v2, min_v2, empReportingC1, ...
    [242,77,155]./255, [0, 139,0]./255, ones(1,3));colormap(cmap_cus2); 
imagesc(simReportingC1); cb2 = colorbar; 
cb2.Ticks = sort([floor(max_v2*100)/100, ceil(min_v2*100)/100, round(empReportingC1,2)]); 
xticks(1:10:length(alpha_sigma));
xticklabels(alpha_sigma(1:10:length(alpha_sigma))); 
yticks(1:10:length(alpha_pC1)); yticklabels(round(alpha_pC1(1:10:length(alpha_pC1)),2));
ylabel('$\alpha_{p_{C=1}}$','interpreter','latex'); 
xlabel('$\alpha_{\sigma}$','interpreter','latex'); %yticks(100:600:2000);
set(gca,'FontSize',15);

%%
param.sigma_deltaT     = 120;
param.sigmaP_deltaT_C1 = 50;
param.sigmaP_deltaT_C2 = 800;
param.mu_P_deltaT_C1   = 0;
param.mu_P_deltaT_C2   = 0;

[propC1_, pCommon_end_] = deal(NaN(length(sigma_A_learning_vec), length(sigma_V_learning_vec),...
    length(alpha_pC1), length(D.numSims)));

for i = 1:length(sigma_A_learning_vec)
    disp(i)
    param.sigma_A_learning = sigma_A_learning_vec(i);
    for j= 1:length(sigma_V_learning_vec)
        param.sigma_V_learning = sigma_V_learning_vec(j);
        for k = 1:length(alpha_pC1)
            for l = 1:D.numSims
                [propC1_(i,j,k,l), pCommon_end_(i,j,k,l), ~] = simulate_learningPhase_wTiming(...
                    param, alpha_pC1(k), D, model);
            end
        end
    end
end









