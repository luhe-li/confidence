%This script aims to run simulations that replicate increases/decreases in
%the common-cause prior after repeated exposure to congruent/incongruent
%audiovisual stimulus pairs
clear all; close all; clc

%% set ground truths
%best-fitting model
model.M_pC1       = 'fullModel'; %4 common-cause priors
model.M_estimates = 'MA';  %model averaging
model.M_unityJdg  = 'MAP'; %estimates-based (heuristic)
%parameters
param.a_A_true    = 2.49;    %proportional bias in auditory spatial perception
param.b_A_true    = -5.60;    %constant bias in auditory spatial perception
param.mu_P        = 0;    %the center of the prior over stimulus location
param.sigma_P     = 100;  %the spread of the prior over stimulus location, 12.49
param.pCommon     = 0.87; %0.27,0.47,0.67, 0.87 common-cause prior one starts with (pre-learning)
param.c_unityJdg  = 3.5;  %the internal criterion 
%the variability of sensory measurements during the learning phase could 
%have a wide range
mean_est_sigma_A  = 9.13;
mean_est_sigma_V  = 1.31;
sigma_A_learning  = linspace( 0.1,  round(mean_est_sigma_A/5*7,2), 21); %logspace(-1, .6, 20); 
sigma_V_learning  = linspace(0.01,  round(mean_est_sigma_V/5*7,2), 21); %logspace(-2, .3, 20);
alpha_pC1         = 0.03;%logspace(-2, -0.6, 5); 
align_error_aA    = param.a_A_true.*[0.75, 0.5];
align_error_bA    = param.b_A_true + [-3, -6]; %0;
align_error       = combvec(align_error_aA, align_error_bA);

%% congruent condition: stimulus pairs presented during the learning phase
sV               = -12:8:12;    %unique visual locations
sA               = (sV - param.b_A_true)./param.a_A_true;
lenS             = length(sV);  %number of visual locations
nT_perLoc        = 40;          %number of trials per location
bool_unity       = [0,0,0,0,1]; %unity judgment was given 20% of trials (1 out of 5)
%pseudo-randomly generate AV pairs
AVpairs_idx      = arrayfun(@(x) randperm(lenS), 1:nT_perLoc, 'UniformOutput', false); 
AVpairs_idx      = cell2mat(AVpairs_idx')';
D.AVpairs(1,:)   = sA(AVpairs_idx(:)); 
D.AVpairs(2,:)   = sV(AVpairs_idx(:));
D.totalTrials    = size(D.AVpairs,2);
D.numSims        = 100; %20 for testing
unity_idx        = arrayfun(@(x) bool_unity(randperm(length(bool_unity))),...
                    1:(D.totalTrials/length(bool_unity)), 'UniformOutput', false);
bool_unity_idx   = cell2mat(unity_idx')';
D.bool_unity     = bool_unity_idx(:)';
%intialize matrices
[propC1, pCommon_end] = deal(NaN(size(align_error_aA,2), length(sigma_A_learning),...
                        length(sigma_V_learning), length(D.numSims)));

%start simulations
for i = 1:size(align_error,2)
    disp(i)
    %unique auditory locations
    param.a_A = align_error(1,i);
    param.b_A = align_error(2,i);
    for j = 1:length(sigma_A_learning)
        disp(j)
        param.sigma_A_learning = sigma_A_learning(j);
        for k = 1:length(sigma_V_learning)
            param.sigma_V_learning = sigma_V_learning(k);
            for l = 1:D.numSims
                [propC1(i,j,k,l), pCommon_end(i,j,k,l), ~] = ...
                    simulate_learningPhase(param, alpha_pC1, D, model);
            end
        end
    end
end

%% plot the figure
delta_pCommon = mean(pCommon_end,4) - param.pCommon;
min_delta     = min(delta_pCommon(:)); 
max_delta     = max(delta_pCommon(:));
cmap_cus      = customize_colormap(max_delta, min_delta, 0); 

figure
for i = 1:size(align_error,2)
    subplot(1,size(align_error,2),i)
    imagesc('XData',sigma_V_learning,'YData',sigma_A_learning,...
        'CData',squeeze(delta_pCommon(i,:,:))); colormap(cmap_cus); 
    caxis([min_delta, max_delta]); hold on
    plot([mean_est_sigma_V,mean_est_sigma_V],[sigma_A_learning(1), ...
        mean_est_sigma_A],'k--','lineWidth',1.5);
    plot([sigma_V_learning(1), mean_est_sigma_V], [mean_est_sigma_A, ...
        mean_est_sigma_A],'k--','lineWidth',1.5); 
    plot(mean_est_sigma_V, mean_est_sigma_A, 'k*', 'lineWidth',1.5,...
        'MarkerSize',10);hold off
    %boundaries, ticks and stuff
    xlim([sigma_V_learning(1), sigma_V_learning(end)]);
    ylim([sigma_A_learning(1), sigma_A_learning(end)]);
    if i == 1
        yticks(round(sort([mean_est_sigma_A, sigma_A_learning([1, 11,...
            length(sigma_A_learning)])]),2)); 
        ylabel('$\sigma_{AV,A}''$ (congruent learning phase)',...
            'interpreter','latex'); 
    else; yticks([]);
    end
    xticks(round(sort([mean_est_sigma_V, sigma_V_learning([1,...
            length(sigma_V_learning)])]),2));
    title(['Measured $a_A = $', num2str(round(align_error(1,i),2)), ...
        '$, b_A = $', num2str(round(align_error(2,i),2))], 'interpreter', 'latex');
    colorbar; caxis([min_delta, max_delta]); 
    xlabel('$\sigma_{AV,V}''$ (congruent learning phase)','interpreter','latex');
    set(gca,'FontSize',15);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.40]);
set(gcf,'PaperUnits','centimeters','PaperSize',[60 15]);
saveas(gcf, ['Simtest_cong_hypSubj_pCommon', num2str(param.pCommon)], 'pdf'); 

%% incongruent condition: stimulus pairs presented during the learning phase
sA_incong     = [sA, sA];
sV_incong     = [(-12:8:12) - 15, (-12:8:12) + 15]; 
lenS_incong   = length(sV_incong);
nT_perLoc     = 20;  %number of trials per location

%pseudo-randomly generate AV pairs
AVpairs_idx           = arrayfun(@(x) randperm(lenS_incong), 1:nT_perLoc, 'UniformOutput', false); 
AVpairs_idx           = cell2mat(AVpairs_idx')';
D_incong.AVpairs(1,:) = sA_incong(AVpairs_idx(:)); 
D_incong.AVpairs(2,:) = sV_incong(AVpairs_idx(:));
D_incong.totalTrials  = size(D_incong.AVpairs,2);
D_incong.numSims      = 100; %20 for testing
D_incong.bool_unity   = D.bool_unity;
%the variability of sensory measurements during the learning phase could 
%have a wide range
sigma_A_learning      = linspace(1, round(mean_est_sigma_A*2,2), 21);%linspace( 1,  8, 20); 
sigma_V_learning      = linspace(0.1, round(mean_est_sigma_V*3,2), 21);%linspace( 0.1,  6, 20); 

%intialize matrices
[propC1_incong, pCommon_end_incong] = deal(NaN(size(align_error_aA,2),...
    length(sigma_A_learning),length(sigma_V_learning), length(D.numSims)));

%start simulations
for i = 1:size(align_error,2)
    disp(i)
    %unique auditory locations
    param.a_A = align_error(1,i);
    param.b_A = align_error(2,i);
    for j = 1:length(sigma_A_learning)
        disp(j)
        param.sigma_A_learning = sigma_A_learning(j);
        for k = 1:length(sigma_V_learning)
            param.sigma_V_learning = sigma_A_learning(k);
            for l = 1:D.numSims
                [propC1_incong(i,j,k,l), pCommon_end_incong(i,j,k,l), ~] = ...
                    simulate_learningPhase(param, alpha_pC1, D_incong, model);
            end
        end
    end
end

%% plot the figure
delta_pCommon_incong = mean(pCommon_end_incong,4) - param.pCommon;
min_delta_incong     = min(delta_pCommon_incong(:)); 
max_delta_incong     = max(delta_pCommon_incong(:));
cmap_cus_incong      = customize_colormap(max_delta_incong, min_delta_incong, 0,...
     [255,165,0]./255, [52,113,131]./255, ones(1,3)); 
% cmap_cus_incong = customize_colormap(max_delta_incong, min_delta_incong, 0); 

figure
for i = 1:size(align_error,2)
    subplot(1,size(align_error,2),i)
    imagesc('XData',sigma_V_learning,'YData',sigma_A_learning,...
        'CData',squeeze(delta_pCommon_incong(i,:,:)));
    colormap(cmap_cus_incong); colorbar; hold on
    plot([mean_est_sigma_V,mean_est_sigma_V],[sigma_A_learning(1), ...
        mean_est_sigma_A],'k--','lineWidth',1.5);
    plot([sigma_V_learning(1), mean_est_sigma_V], [mean_est_sigma_A, ...
        mean_est_sigma_A],'k--','lineWidth',1.5); 
    plot(mean_est_sigma_V, mean_est_sigma_A, 'k*', 'lineWidth',1.5,...
        'MarkerSize',10);hold off
    caxis([min_delta_incong, max_delta_incong]);
    xlim([sigma_V_learning(1), sigma_V_learning(end)]);
    ylim([sigma_A_learning(1), sigma_A_learning(end)]);
    if i == 1
        yticks(round(sort([mean_est_sigma_A, sigma_A_learning([1,...
            length(sigma_A_learning)])]),2)); 
        ylabel('$\sigma_{AV,A}''$ (incongruent learning phase)','interpreter','latex'); 
    else; yticks([]);
    end
    title(['Measured $a_A = $', num2str(round(align_error(1,i),2)), ...
        '$, b_A = $', num2str(round(align_error(2,i),2))], 'interpreter', 'latex');
    xticks(round(sort([mean_est_sigma_V, sigma_V_learning([1,...
            length(sigma_V_learning)])]),2)); 
    xlabel('$\sigma_{AV,V}''$ (incongruent learning phase)','interpreter','latex');
    set(gca,'FontSize',15);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.40]);
set(gcf,'PaperUnits','centimeters','PaperSize',[60 15]);
%saveas(gcf, ['Simtest_incong_hypSubj_pCommon', num2str(param.pCommon)], 'pdf'); 


