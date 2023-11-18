%This script aims to analyze the unity judgments as well as auditory
%localization responses from the bimodal spatial-localization task (AV)
clear all; close all; clc; rng(1)
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'Unimodal localization v3/Data/']));
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'Bimodal localization v2/Data/']));

%subject information
subjNs = [   3,   4,   5,   6,   8,   9,  11,  12,  13,  18,  19,  20,...
                21,  22,  23,  24,  25]; %15AD, 16SM, 17SX are identified as outliers
subjIs = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','ZL','RE','MM',...
              'LL','CS','JH','MD','HHL'};
lenS   = length(subjNs);

%% load data from the unimodal localization task
[Aloc_unique, avgAlocR_uni] = deal(NaN(lenS, 4)); %2 auditory stimulus locations
for s = 1:lenS
    %To compute VE, we need to load unimodal localization responses
    D = load(['Unimodal_localization_sub', num2str(subjNs(s)),'.mat'],...
              'Unimodal_localization_data');
    unimodalA_data   = D.Unimodal_localization_data{end}.data;
    Aloc_uni         = unimodalA_data(1,:);
    AlocR_uni        = unimodalA_data(2,:);
    
    %compute the mean unimodal auditory localization responses
    Aloc_unique(s,:) = unique(Aloc_uni);
    for j = 1:size(Aloc_unique,2)
        idx_A = find(abs(Aloc_uni - Aloc_unique(s,j)) < 1e-1);
        avgAlocR_uni(s,j) = mean(AlocR_uni(idx_A));
    end
end

%% load data from the main experiment
%number of total sessions
numSes        = 2;
%unique visual locations (there are four)
Vloc_unique   = -12:8:12;
%initializations
[nTT_perSes, Aloc, Vloc, AVpairs_order, locR, Ujdg] = deal(cell(1, lenS));

for s = 1:lenS
    nTT_perSes{s} = [0, NaN(1, numSes)];
    %1st dim: in the preceding trial, participants report C=1, C=2
    %2nd dim: in the preceding trial, the spatial discrepancy is small vs. large
    for i = 1:numSes
        C = load(['BimodalLocalization_pre_sub', num2str(subjNs(s)), '_session',...
            num2str(i),'.mat'], 'BimodalLocalization_pre_data');
        %the number of total trials per session
        nTT_perSes{s}(i+1)            = C.BimodalLocalization_pre_data{1}.numTotalTrials;
        %all tested auditory locations (pseudorandomized)
        idx_l                         = sum(nTT_perSes{s}(1:i)) +1;
        idx_r                         = sum(nTT_perSes{s}(1:i+1));
        Aloc{s}(idx_l:idx_r)          = C.BimodalLocalization_pre_data{end-1}.arrangedLocs_deg;
        %all tested visual locations (pseudorandomized)
        Vloc{s}(idx_l:idx_r)          = C.BimodalLocalization_pre_data{end-2}.arrangedLocs_deg;
        %AV pairs were presented in pseudo-random order
        AVpairs_order{s}(idx_l:idx_r) = C.BimodalLocalization_pre_data{1}.AVpairs_order;
        %localizatino responses
        locR{s}(idx_l:idx_r)          = C.BimodalLocalization_pre_data{end}.localization;
        %unity judgments
        Ujdg{s}(idx_l:idx_r)          = C.BimodalLocalization_pre_data{end}.unity;
    end
end

%% splitting the data into different groups and compute separately the 
%ventriloquism effects and the proportion of reporting a common cause 
%ExpInfo.AVpairs_allComb = combvec(1:4, 1:4, 1:2); %1st row: A loc, 2nd row: V loc, 3rd row: localization modality
%      1   2   3   4   1   2   3   4   1   2   3   4   1   2   3   4
%      1   1   1   1   2   2   2   2   3   3   3   3   4   4   4   4
%idx = 1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16

%abs spatial discrepancy = 8 (Group1), 16 (Group2), 24 (Group3)
%            sV  =  -12,      -4,      4,       12
% sA = -12          N/A     Group1   Group2   Group3
% sA = -4          Group1     N/A    Group1   Group2
% sA = 4           Group2   Group1     N/A    Group1
% sA = 12          Group3   Group2   Group1     N/A
%
len_stim  = length(Vloc_unique);
grp_num   = {[2, 5, 7, 10, 12, 15],[3, 8, 9, 14], [4, 13]};
nC        = 2; %number of causes
[locR_cond, propC1_cond] = deal(NaN(lenS, nC, length(Vloc_unique)-1));

for s = 1:lenS
    %splitting trials based on unity-judgment responses (C=1 or C=2)
    idx_predC    = cell(1,nC);
    idx_predC{1} = find(Ujdg{s} == 1);
    idx_predC{2} = find(Ujdg{s} == 2);
    
    %repeat the mean unimodal localization responses so that the
    %vector has the same length as Aloc{s}(idx_cond{k,j}+1)
    avgAlocR_uni_rep_temp = Aloc{s};
    for p = 1:len_stim
        avgAlocR_uni_rep_temp(avgAlocR_uni_rep_temp == Aloc_unique(s,p)) = ...
            avgAlocR_uni(s,p);
    end
    
    %idx_cond saves all the indices for 6 different groups
    %(2 unity judgment responses x 3 absolute spatial discrepancy) 
    idx_cond = cell(nC, length(grp_num));
    for j = 1:length(grp_num)
        %find idx for each group (divided based on spatial discrepancy)
        idx_slc_temp = ismember(AVpairs_order{s}, grp_num{j});
        idx_slc = find(idx_slc_temp);
        
        %preceding trial: C=1 or C=2
        for k = 1:nC 
            idx_cond{k,j} = sort(intersect(idx_slc, idx_predC{k}));
            
            %if idx_cond contains trial 640, get rid of that
            if ismember(idx_r, idx_cond{k,j}); idx_cond{k,j}(end) = [];end
            
            %compute ventriloquism effects
            %idx + 1 for next trial
            avgAlocR_uni_rep = avgAlocR_uni_rep_temp(idx_cond{k,j}+1);   
            locR_cond(s,k,j) = mean((locR{s}(idx_cond{k,j}+1) - avgAlocR_uni_rep).*...
                    -sign(Aloc{s}(idx_cond{k,j}+1)));
            propC1_cond(s,k,j) = sum(Ujdg{s}(idx_cond{k,j}+1) == 1)/...
                length(Ujdg{s}(idx_cond{k,j}+1));
        end
    end
end

%% plot VE
propC1_cond_avg  = squeeze(nanmean(propC1_cond,1));
lenS_notnan      = squeeze(sum(~isnan(propC1_cond)));
propC1_cond_sem  = sqrt(propC1_cond_avg.*(1-propC1_cond_avg)./lenS_notnan);

absSpatialD      = repmat(8:8:24, [2,1]);
lenS_notnan      = squeeze(sum(~isnan(locR_cond)));

locR_cond_avg    = squeeze(nanmean(locR_cond,1));
locR_cond_sem    = squeeze(nanstd(locR_cond,1))./sqrt(lenS_notnan);
loc_cond_avg_pct = locR_cond_avg./absSpatialD.*100;
loc_cond_sem_pct = sqrt((locR_cond_avg./absSpatialD).*(1-locR_cond_avg./absSpatialD)./lenS_notnan).*100;


%variables for plotting
orange          = [255, 126,0]./255;
green           = [76,153,0]./255;
grey            = .5.*ones(1,3);
colorMap1       = [linspace(1,orange(1),255)', linspace(1,orange(2),255)',...
                   linspace(1,orange(3),255)'];
colorMap2       = [linspace(1,green(1),255)', linspace(1,green(2),255)',...
                   linspace(1,green(3),255)'];
colorMap3       = [linspace(1,grey(1),255)', linspace(1,grey(2),255)',...
                   linspace(1,grey(3),255)'];             

figure;
subplot(1,3,1)
h(1) = heatmap(locR_cond_avg,'XDisplayLabels', {'8','16','24'}, 'YDisplayLabels', {'C=1','C=2'});
%xlabel('Absolute spatial discrepancy (deg) in the preceding trial');
%ylabel('Unity judgment in the preceding trial');
colormap(h(1),colorMap1); colorbar; caxis([floor(min(locR_cond_avg(:))), ...
    ceil(max(locR_cond_avg(:)))]); grid off; set(gca,'FontSize',15);
ylabel(sprintf('Unity judgment \nin the preceding trial'));
title('Ventriloquims effects (deg)'); set(gca,'FontSize',15);

subplot(1,3,2)
h(2) = heatmap(loc_cond_avg_pct,'XDisplayLabels', {'8','16','24'}, 'YDisplayLabels', {'C=1','C=2'});
xlabel('Absolute spatial discrepancy (deg) in the preceding trial');
%ylabel('Unity judgment in the preceding trial');
colormap(h(2), colorMap2); colorbar; caxis([floor(min(loc_cond_avg_pct(:))), ...
    ceil(max(loc_cond_avg_pct(:)))]); grid off;
title(sprintf('Ventriloquims effects divided by \nabsolute spatial discrepancy (%%)')); 
set(gca,'FontSize',15);

subplot(1,3,3)
h(3) = heatmap(propC1_cond_avg,'XDisplayLabels', {'8','16','24'}, 'YDisplayLabels', {'C=1','C=2'});
%xlabel('Absolute spatial discrepancy (deg) in the preceding trial');
%ylabel('Unity judgment in the preceding trial');
colormap(h(3), colorMap3); colorbar; caxis([0,1]); grid off; 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.9,0.3]);
set(gcf,'PaperUnits','centimeters','PaperSize',[55 12]);
title('The probability of reporting C=1'); set(gca,'FontSize',15);
%saveas(gcf, 'VE_propC1_cond_heatmap_allSubjs', 'pdf'); 

%% add error bars
jitter = [-0.05, 0.05];

figure
subplot(1,3,1)
for i = 1:nC
    errorbar((1:3) + jitter(i), locR_cond_avg(i,:), locR_cond_sem(i,:),'lineWidth',2); hold on
end
xlim([0.5,3.5]); ylim([2,7]); box off; grid on
xticks(1:3); xticklabels({'8','16','24'}); legend({'C=1','C=2'});
ylabel('Ventriloquims effects (deg)');
set(gca,'FontSize',15);

subplot(1,3,2)
for i = 1:nC
    errorbar((1:3) + jitter(i), loc_cond_avg_pct(i,:), loc_cond_sem_pct(i,:),'lineWidth',2); hold on
end
xlim([0.5,3.5]); ylim([0, 70]); box off; grid on
xticks(1:3); xticklabels({'8','16','24'}); 
xlabel(sprintf('Absolute spatial discrepancy (deg) in the preceding trial'));
ylabel(sprintf('Ventriloquim effects \ndivided by absolute \nspatial discrepancy (%%)'));
set(gca,'FontSize',15);

subplot(1,3,3)
for i = 1:nC
    errorbar((1:3) + jitter(i), propC1_cond_avg(i,:), propC1_cond_sem(i,:),'lineWidth',2); hold on
end
xlim([0.5,3.5]); ylim([0,1]); box off; grid on
xticks(1:3); xticklabels({'8','16','24'}); 
%xlabel(sprintf('Absolute spatial discrepancy (deg) in the preceding trial'));
ylabel(sprintf('The probability of \nreporting C=1'));
set(gca,'FontSize',15);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7,0.3]);
set(gcf,'PaperUnits','centimeters','PaperSize',[55 12]);
saveas(gcf, 'VE_propC1_cond_errorbar_allSubjs', 'pdf'); 

