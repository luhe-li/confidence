%This script aims to compare the standard deviation of localization 
%responses from the unimodal spatial localization task as a function of 
%stimulus location
clear all; close all; clc; rng(1)
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'Unimodal localization v3/Data']));
subjNs = [   3,   4,   5,   6,   8,   9,  11,  12,  13, ...
            18,  19,  20,  21,  22,  23,  24,  25];
subjIs = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR',...
          'ZL','RE','MM','LL','CS','JH','MD','HHL'};
lenS   = length(subjNs);
Vloc   = -12:8:12;

%% 
[var_locError_A, var_locError_V] = deal(NaN(lenS, 4));
for i = 1:lenS
    subjN    = subjNs(i);
    subjI    = subjIs{i}; 
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'Unimodal localization v3/Data/',subjI]));
    C        = load(['Unimodal_localization_sub', num2str(subjN), '.mat'],...
                'Unimodal_localization_data');
    Adata     = C.Unimodal_localization_data{end}.data;
    sA_unique = unique(Adata(1,:));
    Vdata     = C.Unimodal_localization_data{end-1}.data;
    for j = 1:length(Vloc)
        idx_A = find(Adata(1,:) == sA_unique(j));
        locResp_A = Adata(2,idx_A);
        var_locError_A(i,j) = std(locResp_A);
        idx_V = find(Vdata(1,:) == Vloc(j));
        locResp_V = Vdata(3,idx_V);
        var_locError_V(i,j) = std(locResp_V);
    end
end

%% plot it
mean_varA_allSubjs = mean(var_locError_A);
sem_varA_allSubjs  = std(var_locError_A)./sqrt(lenS);
mean_varV_allSubjs = mean(var_locError_V);
sem_varV_allSubjs  = std(var_locError_V)./sqrt(lenS);

figure
errorbar(Vloc, mean_varA_allSubjs, sem_varA_allSubjs,'b','lineWidth',2); hold on
errorbar(Vloc, mean_varV_allSubjs, sem_varV_allSubjs,'r','lineWidth',2); hold off; box off
xlabel('Stimulus location (dva)'); ylabel('SD of localization errors');
xticks(Vloc); set(gca,'FontSize',15);






