clear all; close all; clc
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'ModelFitting/Figures_results']));
addpath(genpath('/Users/hff/Desktop/NYU/Project2/Experiment code/ModelFitting/Fits'));
%subject number (outlier participants: 15, 16, 17)
subjN_dict      = [   3,   4,   5,   6,   8,   9,  11,  12,  13,...
                     18,  19,  20,  21,  22,  23,  24, 25];
%subject initial (outlier participants: 'AD', 'SM', 'SX'
subjI_dict      = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR',...
                   'ZL','RE','MM','LL','CS','JH','MD','HHL'};      
lenS            = length(subjN_dict);
             
%% load the file
idx_mat    = {[9,10;9,11],[9,10;11,12],[9,9;9,9]};
C = load('Results_modelComparison_overall.mat', 'ModelComparison');
bestM = C.ModelComparison{end};
[sigma_AV_A, sigma_AV_V] = deal(NaN(1, lenS));
pC1_allConds = NaN(lenS, 2,2); %2 conditions; 2 phases
delta_pC1_allConds = NaN(lenS, 2);
for s = 1:lenS
    bestM_s   = bestM{s};
    seg_idx   = find(bestM_s == '-'); 
    strategy1 = bestM_s((seg_idx(1)+1):(seg_idx(2)-1));
    strategy2 = bestM_s((seg_idx(2)+1):end);    
    D         = load(['ModelFitting_updatePrior_', strategy1,'_',strategy2,...
                    '_sub', num2str(subjN_dict(s)), '.mat'], 'ModelFitting');
    sigma_AV_A(s) = D.ModelFitting{end}.P_f(6);
    sigma_AV_V(s) = D.ModelFitting{end}.P_f(7);
    
    switch strategy1
        case 'HighPlasticityLongLasting'; idx_mat_s = idx_mat{2};
        case 'HighPlasticityShortLasting'; idx_mat_s = idx_mat{1};
        case 'NoPlasticity'; idx_mat_s = idx_mat{3};
    end
    
    for i = 1:2 %2 conditions
        pC1_allConds(s,i,:) = D.ModelFitting{end}.P_f(idx_mat_s(i,:));
        delta_pC1_allConds(s,i) = pC1_allConds(s,i,2) - pC1_allConds(s,i,1);
    end
end

%%
figure
scatter(sigma_AV_A, delta_pC1_allConds(:,1)',180,'filled', 'r','Marker',...
    'o','MarkerFaceAlpha',0.2,'MarkerEdgeColor','r'); hold on;
scatter(sigma_AV_A, delta_pC1_allConds(:,2)',180,'filled', 'b','Marker',...
    'square','MarkerFaceAlpha',0.2,'MarkerEdgeColor','b'); hold off; box off; grid on
xlim([0,30]); ylim([-0.4,0.6]); 
xlabel(sprintf('The variability of sensory measurements \nduring the pre- & post-learning phases'));
ylabel('Change of the common-cause prior');legend({'Congruent','Incongruent'});
set(gca,'FontSize',15);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.3, 0.5]);
set(gcf,'PaperUnits','centimeters','PaperSize',[25 25]);
%saveas(gcf,'Delta_pC1_vs_sigma_AV_A.pdf');



