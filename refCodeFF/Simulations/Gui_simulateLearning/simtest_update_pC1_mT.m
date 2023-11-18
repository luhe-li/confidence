%% set ground truth
clear all; close all; clc
%load parameter estimate
subjN = 6; subjI = 'YZ'; sesNum = 2; %incongruent condition
addpath(genpath(['/Users/Guinevere/Desktop/NYU/Project2/ModelFitting/Fits/',subjI]));
C = load(['ModelFitting_updatePrior_fullModel_unityJdg1_locResp1_sub', num2str(subjN), '.mat'],...
    'ModelFitting');
addpath(genpath(['/Users/Guinevere/Desktop/NYU/Project2/Adaptation v2/Data/',subjI]));
F = load(['Adaptation_incongruent_sub', num2str(subjN),'_session',num2str(sesNum),'.mat'],...
    'Adaptation_incongruent_data');
estimatedP                = C.ModelFitting{end}.P_f;
param.sigma_DeltaT        = 120; %120
param.mu_P_DeltaT_C1      = 0;
param.sigma_P_DeltaT_C1   = 50;
param.mu_P_DeltaT_C2      = 0; %250ms if assume a peripheral prior
param.sigma_P_DeltaT_C2   = 160; %linspace(100,2000,10); FREE PARAMETER
param.sigma_P_spatial     = 100;
param.mu_P_spatial        = 0;
param.criterion           = 0.5;
param.sigma_spatial_AV_A  = 10;%estimatedP(5); %deg
param.sigma_spatial_AV_V  = estimatedP(6);
param.a_A                 = estimatedP(1);
param.b_A                 = estimatedP(2);
param.pCommon             = estimatedP(9);
alpha                     = linspace(1e-4, 0.3, 10);
alpha_deltaT              = linspace(1e-2, 0.3, 10);

D.AVpairs(1,:)        = F.Adaptation_incongruent_data{4}.arrangedLocs_deg; %A loc
D.AVpairs(3,:)        = F.Adaptation_incongruent_data{3}.arrangedLocs_deg; %V loc
D.AVpairs(2,:)        = -F.Adaptation_incongruent_data{3}.timing_relative./2; %A time relative
D.AVpairs(4,:)        = F.Adaptation_incongruent_data{3}.timing_relative./2; %V time relative
D.totalTrials         = size(D.AVpairs,2);
D.numSims             = 1e2;
[propC1, pCommon_end] = deal(NaN(length(alpha_deltaT),length(alpha),length(D.numSims)));

for i = 1:length(alpha_deltaT)
    disp(i)
    param.alpha_deltaT = alpha_deltaT(i);
    for j = 1:length(alpha)
        param.alpha = alpha(j);
        for k = 1:D.numSims
            [propC1(i,j,k), pCommon_end(i,j,k)] = simulateDissociationPhase_update_pC1_mT(param, D);
        end
    end
end

%%
simChange_pCommon = mean(pCommon_end,3)-param.pCommon;
empChange_pCommon = estimatedP(10) - param.pCommon;
%tol = 0.05;
%bool_withinTol = abs(simChange_pCommon-empChange_pCommon) < tol;
%simChange_pCommon(~bool_withinTol) = 0;
simReportingC1 = mean(propC1,3);
empReportingC1 = sum(F.Adaptation_incongruent_data{end}.unity==1)/...
    sum(~isnan(F.Adaptation_incongruent_data{end}.unity));
% tol = 0.1;
% bool_withinTol = abs(simReportingC1-empReportingC1) < tol;
% simReportingC1(~bool_withinTol) = 0;
figure(1)
imagesc(alpha, alpha_deltaT, simChange_pCommon); colorbar
xlabel('\alpha_{p_{C=1}}'); ylabel('\alpha_{\Delta_t}'); %yticks(100:600:2000);
set(gca,'FontSize',15);
    
figure(2)
imagesc(alpha, alpha_deltaT, simReportingC1); colorbar
xlabel('\alpha_{p_{C=1}}'); ylabel('\alpha_{\Delta_t}'); %yticks(100:600:2000);
set(gca,'FontSize',15);



