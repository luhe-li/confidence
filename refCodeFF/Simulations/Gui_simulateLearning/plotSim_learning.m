%% plot simulations with estimated pCommone
clear all; close all; clc
subjNs     = [3,4,5,6,8,9,11,12,13,15,16,17,18];
subjIs     = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','AD','SM','SX','ZL'};
cond_order = [2,1;1,2;2,1;1,2;2,1;1,2;2,1;2,1;1,2;2,1;1,2;1,2;2,1];%1: cong; 2: incong

%%
global param_name grid_name Cond
Cond       = {'congruent','incongruent'};
param_name = {'sigma_deltaT','muP_deltaT_C1','sigmaP_deltaT_C1',...
                'muP_deltaT_C2','sigmaP_deltaT_C2','alpha_pC1_congruent',...
                'alpha_pC1_incongruent'};
grid_name  = {'sigma_AV_A_lb_ub_congruent','sigma_AV_V_lb_ub_congruent',...
                'sigma_AV_A_lb_ub_incongruent','sigma_AV_V_lb_ub_incongruent'};
% cMAP  = [125,203,151;89,191,182;78,135,163;79,79,119;93,46,88;...
%         148,48,81;195,56,79;243,115,83;245,159,77;249,204,84;237,225,121;...
%         187,216,123;210,105,30]./255;
cMAP =  [242,182,126;100,195,187;215,133,76;91,144,169]./255;

%% load simulations
subjN_selected = [16,6,15,8];%[6,8,9,11];
cond_selected  = [1,1,2,2];%[1,2,1,2];
subjI_selected = {'SM','YZ','AD','NH'};%{'YZ','NH','ZZ','BB'};
bool_save      = 1;
for i = 1:length(subjN_selected)
    idx_subj   = find(subjN_selected(i) == subjNs);
    sesNum     = cond_order(idx_subj,:); %incongruent condition
    plt_simPC1(subjN_selected(i), subjI_selected{i}, idx_subj, ...
        sesNum,cond_selected(i), cMAP(i,:), bool_save);
end


