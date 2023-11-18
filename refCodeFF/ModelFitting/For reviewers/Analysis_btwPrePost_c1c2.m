%% get needed info
clear all; close all; clc
%subject number (outlier participants: 15, 16, 17)
subjN_dict = [   3,   4,   5,   6,   8,   9,  11,  12,  13,  15,  16,  17,...
                18,  19,  20,  21,  22,  23,  24, 25];
%subject initial (outlier participants: 'AD', 'SM', 'SX'
subjI_dict = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','AD','SM','SX',...
              'ZL','RE','MM','LL','CS','JH','MD','HHL'};
%the order of condition was counterbalanced across participants
%1: congruent; 2: incongruent
order_dict = [ 2,1; 1,2; 2,1; 1,2; 2,1; 1,2; 2,1; 2,1; 1,2; 2,1; 1,2; 1,2;...
               2,1; 1,2; 2,1; 2,1; 1,2; 2,1; 1,2; 2,1]; 
%create popout window
prompt     = {'Subject ID:','Plot unity judgment (1:yes; 0: no):',...
                'Plot demeand locR:', 'Plot mean locR:',...
                'save data', 'export tables for percent C1 and VE',...
                'export a table for unityJdg', 'export a table for locR:', ...
                'export a table for locR (abs(spatialD)):'};
dlgtitle   = 'Input';
dims       = [1 35];
definput   = {'3','1','1','1','0','0','0','0','0'};
answer     = inputdlg(prompt,dlgtitle,dims,definput);
bool_plt   = arrayfun(@(idx) str2double(answer(idx)), 2:4); 
bool_save  = arrayfun(@(idx) str2double(answer(idx)), 5:length(prompt)); 

%% load data
%which participant is selected
subjN = str2double(answer(1));
subjI = subjI_dict{find(subjN==subjN_dict,1)};
%load data
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                    'Unimodal localization v3/Data/', subjI]));
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                    'Bimodal localization v2/Data/', subjI]));
cond     = {'congruent','incongruent'};
order    = order_dict(find(subjN==subjN_dict,1),:);
            %[1,2]: congruent first and incongruent second
            %[2,1]: incongruent first and congruent second
%number of condition: congruent, incongruent
lenC     = length(cond);
%phase of testing the spatial perception
phase    = {'pre','post'}; 
%number of phases
lenP     = length(phase);
%response modality
modality = {'A','V'};
%number of response modality
lenM     = length(modality);
%number of total trials for each phase
%4 x 4 = 16 different AV pairs, each pair was repeated 20 times
nTT      = 320; 

%initialize matrices
%AVpairs_order   : shuffled 1-32 with each trial type repeated 10 times
%unityResp       : reported common-cause judgment (1: C=1, 2: C=2)
%localizeModality: cued modality after stimulus presentation (1: A, 2, V)
%locResp         : localization responses (in deg)
[AVpairs_order, unityResp, localizeModality, locResp, unityResp_c1,...
    localizeModality_c1, locResp_c1] = deal(NaN(lenC,lenP, nTT));

%load files and get data we need
for i = 1:lenC %for each condition
    for j = 1:lenP %for each phase
        C = load(['BimodalLocalization_', phase{j},'_sub', num2str(subjN),...
            '_session',num2str(order(i)),'.mat'],['BimodalLocalization_',...
            phase{j}, '_data']);
        AVpairs_order(i,j,:)    = eval(['C.BimodalLocalization_',phase{j},...
                                    '_data{1}.AVpairs_order']); 
        unityResp(i,j,:)        = eval(['C.BimodalLocalization_',phase{j},...
                                    '_data{end}.unity']); 
        localizeModality(i,j,:) = eval(['C.BimodalLocalization_',phase{j},...
                                    '_data{1}.localize_modality']); 
        locResp(i,j,:)          = eval(['C.BimodalLocalization_',phase{j},...
                                    '_data{end}.localization']);
        
        %divide data based on whether the preceding trials is c=1 or c=2                        
        idx_c1                      = find(unityResp(i,j,:) == 1)';
        if ismember(nTT, idx_c1); idx_c1(end) = []; end %exclude the last trial
        idx_next_c1                 = idx_c1 + 1;
        len_idx                     = length(idx_c1);
        unityResp_c1(i,j,1:len_idx) = squeeze(unityResp(i,j,idx_c1))';
        localizeModality_c1(i,j,1:len_idx) = squeeze(localizeModality(i,j,idx_c1))';
        locResp_c1(i,j,1:len_idx)          = squeeze(locResp(i,j,idx_c1))';
    end
end
%other useful info
%different AV stimulus pairs
%1st row: A locs; 2nd row: V locs; 3rd row: response modality
AVpairs_allComb    = C.BimodalLocalization_post_data{1}.AVpairs_allComb;
%number of trial types: 4 (A locs) x 4 (V locs) x 2 (response modality) = 32
numDiffAVpairs     = size(AVpairs_allComb,lenM);
%number of trials per trial type = 10
nT_perPair         = nTT/numDiffAVpairs; 
%auditory location: participant-specific in physical space (dva)
A_loc              = C.BimodalLocalization_post_data{4}.Distance; 
%visual location: -12,-4,4,12
V_loc              = C.BimodalLocalization_post_data{3}.Distance; 
%combination of AV stimulus pairs
%1st row: A loc in perceptual space (-12,-4,4,12)
%2nd row: V loc in perceptual/physical space (-12,-4,4,12)
AV_loc_comb        = combvec(V_loc, V_loc); 
%V-A = [-24, -16, -8, 0, 8, 16, 24] in the scale of visual-equivalent location
AV_discrepancy     = sort(unique(AV_loc_comb(2,:) - AV_loc_comb(1,:))); 
%number of audiovisual spatial discrepancy = 7
lenD               = length(AV_discrepancy); 
%absolute spatial discrepancy
abs_AV_discrepancy = unique(abs(AV_discrepancy));
%number of absolute spatial discrepancy
lenD_abs           = length(abs_AV_discrepancy);

%% load data from unimodal spatial-localization task
C              = load(['Unimodal_localization_sub', num2str(subjN), '.mat'],...
                    'Unimodal_localization_data');
%get the localization responses 
%1st row: physical stimulus location
%2nd row: localization responses (120 trials for each modality)
%3rd row: response time
LocR_A_uni     = C.Unimodal_localization_data{end}.data;
LocR_V_uni     = C.Unimodal_localization_data{end-1}.data;
%compute the mean unimodal localization responses
meanLocR_A_uni = arrayfun(@(idx) mean(LocR_A_uni(2,abs(LocR_A_uni(1,:) - ...
                    A_loc(idx)) < 1e-3)), 1:length(A_loc));
meanLocR_V_uni = arrayfun(@(idx) mean(LocR_V_uni(3,abs(LocR_V_uni(1,:) - ...
                    V_loc(idx)) < 1e-3)), 1:length(V_loc));
meanLocR_uni   = {meanLocR_A_uni, meanLocR_V_uni};

%% organize data
%--------------------------initialize matrices-----------------------------
%unityMat  : stores common-cause judgments
%            2 (conditions) x 2 (phases) x 4 (A locations) x 4 (V locations) x 20 (trials)
%pC1_resp  : stores the percentage of C = 1 responses for each AV pair
%            2 (conditions) x 2 (phases) x 4 (A locations) x 4 (V locations)
%locRespMat: stores the localization responses
%            2 (conditions) x 2 (phases) x 4 (A locations) x 4 (V locations)
%                           x 2 (localization modality) x 10 (trials)
unityMat   = NaN(lenC, lenP, length(A_loc), length(V_loc), nT_perPair*lenM);
pC1_resp   = NaN(lenC, lenP, length(A_loc), length(V_loc));
locRespMat = NaN(lenC, lenP, length(A_loc), length(V_loc), lenM, nT_perPair); 

%pC1_resp_lenD: stores the percentage of C = 1 responses for each spatial discrepancy
%mean_pC1_lenD: nanmean(pC1_resp_lenD)
%-------------------------Spatial discrepancy(V-A)-----------------------------
%       -24,       -16,       -8,        0,         8,        16,        24
%------------------------------------------------------------------------------
%V1: p(_|A4,V1) p(_|A3,V1)  p(_|A2,V1)  p(_|A1,V1)     NaN       NaN       NaN
%V2:    NaN     p(_|A4,V2)  p(_|A3,V2)  p(_|A2,V2)  p(_|A1,V2)   NaN       NaN  
%V3:    NaN         NaN     p(_|A4,V3)  p(_|A3,V3)  p(_|A2,V3) p(_|A1,V3)  NaN   
%V4:    NaN         NaN        NaN      p(_|A4,V4)  p(_|A3,V4) p(_|A2,V4) p(_|A1,V4)
%
%mean_LocResp: stores the mean localization responses (m=1: A)
%V1: avg(r_A4)  avg(r_A3)  avg(r_A2)   avg(r_A1)       0          0        0
%V2:     0      avg(r_A4)  avg(r_A3)   avg(r_A2)   avg(r_A1)      0        0
%V3:     0        0        avg(r_A4)   avg(r_A3)   avg(r_A2)  avg(r_A1)    0
%V4:     0        0            0       avg(r_A4)   avg(r_A3)  avg(r_A2)  avg(r_A1)
%(m=2:V)
%A1: avg(r_V4)  avg(r_V3)  avg(r_V2)   avg(r_V1)       0          0        0
%A2:     0      avg(r_V4)  avg(r_V3)   avg(r_V2)   avg(r_V1)      0        0
%A3:     0          0      avg(r_V4)   avg(r_V3)   avg(r_V2)  avg(r_V1)    0
%A4:     0          0          0       avg(r_V4)   avg(r_V3)  avg(r_V2) avg(r_V1)

pC1_resp_lenD                = NaN(lenC, lenP, length(V_loc), lenD); 
[mean_pC1_lenD, SD_pC1_lenD] = deal(NaN(lenC, lenP, lenD));
[mean_LocResp,  SD_LocResp]  = deal(zeros(lenC, lenP, lenM,length(V_loc), lenD));
%define a function that calculates the SD for bernoulli distribution
calcSD                       = @(p,n) sqrt(p.*(1-p)./n);
  
%loop through all conditions and all phases
for i = 1:lenC %for each condition
    for j = 1:lenP %for each phase
        for k = 1:size(AV_loc_comb,2) %16 AV pairs
            %first find the indices that correspond to the same AV pair
            %e.g., A = V = -12 deg, k2 = [1 (localize A), 17 (localize V)];
            k2 = [k, k + size(AV_loc_comb,2)];
            %find the index that indicates the type of trial
            idx_samePair_localizeA = (AVpairs_order(i,j,:) == k2(1));
            idx_samePair_localizeV = (AVpairs_order(i,j,:) == k2(2));
            unityResp_samePair     = [squeeze(unityResp(i,j,idx_samePair_localizeA));...
                                      squeeze(unityResp(i,j,idx_samePair_localizeV))];
            unityMat(i,j,AVpairs_allComb(1,k), AVpairs_allComb(2,k),:) = ...
                unityResp_samePair;
            pC1_resp(i,j,AVpairs_allComb(1,k), AVpairs_allComb(2,k)) = ...
                sum(unityResp_samePair==1)/(nT_perPair*2);
            %store the localization responses separately for each
            %localization modality
            locRespMat(i,j,AVpairs_allComb(1,k), AVpairs_allComb(2,k),1,:) = ...
                squeeze(locResp(i,j,idx_samePair_localizeA));
            locRespMat(i,j,AVpairs_allComb(1,k), AVpairs_allComb(2,k),2,:) = ...
                squeeze(locResp(i,j,idx_samePair_localizeV));
        end
        
        %mean p(reporting 'C=1') as a function of spatial discrepancy
        for l = 1:length(V_loc) %for each visual/auditory stimulus location
            %l=1: idx=1:4; l=2: idx=2:5; l=3: idx=3:6; l=4: idx=4:7;
            idx_lb = l; idx_ub = floor(lenD/2)+l;
            %need to use fliplr.m because the order is 
            %p('C=1'|A4,V1), p('C=1'|A3,V1), p('C=1'|A2,V1), p('C=1'|A1,V1)
            pC1_resp_lenD(i,j,l,idx_lb:idx_ub) = fliplr(squeeze(pC1_resp(i,j,:,l))');
            for m = 1:lenM %for each modality
                if m == 1; locRespMat_ijm = squeeze(locRespMat(i,j,:,l,m,:))';
                else; locRespMat_ijm = squeeze(locRespMat(i,j,l,:,m,:))'; end
                mean_LocResp(i,j,m,l,idx_lb:idx_ub) = ...
                        squeeze(mean_LocResp(i,j,m,l,idx_lb:idx_ub))'+...
                        fliplr(sum(locRespMat_ijm./nT_perPair));
                SD_LocResp(i,j,m,l,idx_lb:idx_ub) = std(locRespMat_ijm);
            end
        end
        
        %calculate the mean percentage of reporting C=1 as a function of discrepancy
        mean_pC1_lenD(i,j,:) = squeeze(nanmean(pC1_resp_lenD(i,j,:,:)));
        if i==1 && j==1 %numT_AV = [20,40,60,80,60,40,20]; trials for each discrepancy
            numT_AV = nT_perPair*lenM.*sum(~isnan(squeeze(pC1_resp_lenD(i,j,:,:)))); 
            %number of trials for each level of absolute spatial discrepancy (80,120,80,40)
            numT_abs_AV = numT_AV(ceil(lenD/2):end) + [0,fliplr(numT_AV(1:floor(lenD/2)))]; 
        end
        SD_pC1_lenD(i,j,:) = calcSD(squeeze(mean_pC1_lenD(i,j,:))', numT_AV); 
        %display the proportion of reporting common-cause
        fprintf(['Cond: ', cond{i}, ', phase: ', phase{j}, '\n']); 
        disp(squeeze(mean_pC1_lenD(i,j,:))');
    end
end

%% plot the unity judgment as a function of spatial discrepancy
x_bds          = [AV_discrepancy(1) - 2, AV_discrepancy(end) + 2]; 
y_bds          = [-0.05,1.05];
x_ticks_AVdiff = AV_discrepancy; 
y_ticks        = 0:0.25:1;
lw             = 2.5; %lineWidth
fs_lbls        = 25;
fs_lgds        = 20;
jitter_prepost = [1, -1];
cMap_unity     = [0.65, 0.65,0.65;0.1,0.1,0.1];
lgd_cond_phase = {'Pre-learning', 'Post-learning'};
%lgd_pos        = [0.35 0.15 0.05 0.1; 0.80 0.15 0.05 0.1]; 

if bool_plt(1) == 1
    figure(1)
    for i = 1:lenC
        subplot(1,2,i)
        addBackground(x_bds, y_bds, x_ticks_AVdiff, y_ticks)
        for j = 1:lenP
            for m = 1:lenD
                plt(j) = errorbar(AV_discrepancy(m)+jitter_prepost(j), ...
                    mean_pC1_lenD(i,j,m), SD_pC1_lenD(i,j,m),'-o',...
                    'MarkerSize',sqrt(numT_AV(m)).*2,'Color',cMap_unity(j,:),...
                    'MarkerFaceColor',cMap_unity(j,:),'MarkerEdgeColor',...
                    cMap_unity(j,:),'lineWidth',lw,'CapSize',15); hold on
            end
            plot(AV_discrepancy+jitter_prepost(j), squeeze(mean_pC1_lenD(i,j,:)),...
                '-','lineWidth',lw,'Color',cMap_unity(j,:)); hold on;
        end
        %add legends
        if i == 1;text(x_bds(1) + 0.2, 1.02, subjI,'FontSize',fs_lgds); hold on;end
        text(x_bds(1) + 5, 0.01, ['Condition: ', cond{i}], 'FontSize',fs_lgds); hold on
        legend([plt(1) plt(2)], lgd_cond_phase,'Location', 'northeast',...
            'FontSize',fs_lgds); legend boxoff;
        xticks(x_ticks_AVdiff); xlim(x_bds); xlabel('Spatial discrepancy (V - A, deg)'); 
        ylabel(sprintf(['The probability \n of reporting a common cause'])); 
        yticks(y_ticks);ylim(y_bds);
        set(gca,'FontSize',fs_lbls);
    end
    %best size
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.75, 0.60]);
    set(gcf,'PaperUnits','centimeters','PaperSize',[50 25]);
    %saveas(gcf, ['UnityJdg_btwPrePost_',subjI], 'pdf'); 

    %for manuscript
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.40]);
    %set(gcf,'PaperUnits','centimeters','PaperSize',[40 20]);
    %saveas(gcf, ['UnityJdg_btwPrePost_',subjI,'_'], 'pdf'); 
end

%% localization responses as a function of stimulus location 
idx_loc_lb      = arrayfun(@(idx) find(~isnan(pC1_resp_lenD(1,1,:,idx)),1),1:lenD);
                                                      %[4,3,2,1,1,1,1];
idx_loc_ub      = arrayfun(@(idx) ceil(lenD/2)-find(~isnan(...
                      flipud(squeeze(pC1_resp_lenD(1,1,:,idx)))),1), 1:lenD)+1;
                                                      %[4,4,4,4,3,2,1];             
idx_sloc_lb     = fliplr(idx_loc_lb(1,:));            %[1,1,1,1,2,3,4];
idx_sloc_ub     = fliplr(idx_loc_ub(1,:));            %[1,2,3,4,4,4,4];
ybd_locR        = [min(mean_LocResp(:))-max(SD_LocResp(:)),...
                    max(mean_LocResp(:))+max(SD_LocResp(:))]; 
y_ticks_locR    = V_loc;
lgd_label       = cell(1,lenD);for l=1:lenD;lgd_label{l}=num2str(AV_discrepancy(l));end
x_lbl           = {'Auditory','Visual'};
cMap_locR(1,:,:)= [0.2706, 0.3059, 0.5529;0.3216, 0.4824, 0.6471;0.3725, 0.6549, 0.7412;...
                   0.0588, 0.7373, 0.7451;0.3529, 0.7373, 0.4588;0.5686, 0.8039, 0;...
                   0.7373, 0.8471, 0];
cMap_locR(2,:,:)= [0.4275, 0.1059, 0.1647; 0.6157, 0.1490, 0.3843; 0.7216, 0.4078, 0.4824;...
                   0.8667, 0.5686, 0.5725; 0.9451, 0.6353, 0.5333; 0.9490, 0.7098, 0.5020;...
                   0.9608, 0.7922, 0.4431];
subplot_idx     = [1,3;2,4];
                
if bool_plt(2) == 1
    for m = 1:lenM
        figure(5+m)
        %x boundaries change based on modality
        xbd     = [eval([modality{m}, '_loc(1)'])-3, eval([modality{m},'_loc(end)'])+3]; 
        x_ticks = round(eval([modality{m}, '_loc']),1);
        for i = 1:lenC
            for j = 1:lenP
                subplot(lenC, lenP, subplot_idx(i,j))
                addBackground(xbd, ybd_locR, [xbd(1),x_ticks, xbd(end)], ...
                    [ybd_locR(1),y_ticks_locR,ybd_locR(end)]);
                plot(eval([modality{m},'_loc']), V_loc,'Color',ones(1,3).*0.5,...
                    'lineWidth',lw,'lineStyle','--'); hold on
                for k = 1:lenD
                    Eb(k) = errorbar(eval([modality{m},'_loc(idx_sloc_lb(k):idx_sloc_ub(k))']), ...
                        squeeze(mean_LocResp(i,j,m,idx_loc_lb(k):idx_loc_ub(k),k)),...
                        squeeze(SD_LocResp(i,j,m,idx_loc_lb(k):idx_loc_ub(k),k)),'-s',...
                        'Color',cMap_locR(m,k,:),'lineWidth', lw,'MarkerSize',15,...
                        'MarkerEdgeColor',cMap_locR(m,k,:),'MarkerFaceColor',...
                        cMap_locR(m,k,:)); hold on
                end
                hold off; box off; axis square;
                xticks(x_ticks); xlim(xbd); xlabel([x_lbl{m},' stimulus location (deg)']);
                yticks(y_ticks_locR); ylim(ybd_locR); 
                ylabel(sprintf([x_lbl{m},' mean localization\nresponses (deg)']));
                if i==1 && j == 2
                    lgd = legend(Eb,lgd_label,'Location','northeast','FontSize',fs_lgds,...
                        'Orientation','horizontal');legend boxoff; htitle = get(lgd,'Title'); 
                    set(htitle,'String',['Spatial discrepancy (V-A, deg)']);
                end
                if i==1 && j==1; text(xbd(1)+0.2,ybd_locR(end)-1,subjI,'FontSize',fs_lgds); end
                %title(['Condition: ', cond{i}, ', phase: ', phase{j}],'FontSize',fs_lgds);
                set(gca,'FontSize',fs_lbls);
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.85]);
            end
        end
        set(gcf,'PaperUnits','centimeters','PaperSize',[40 40]);
        %saveas(gcf, ['mean_LocResp',modality{m},'_',cond{i}, '_',subjI], 'pdf'); 
    end
end

%% compute ventriloquism effects
%initialize matrices
meanLocR_uni_rep{1}    = repmat(meanLocR_uni{1},[length(V_loc), 1]);
meanLocR_uni_rep{2}    = repmat(meanLocR_uni{2}',[1, length(A_loc)]);
[meanVE,sdLocR_pooled,seLocR_pooled] = deal(NaN(lenC, lenP, lenM,lenD)); %positive means shifting to V loc
locR                   = cell(lenC, lenP, lenM, lenD);
%we combine data that have the same absolute spatial discrepancy
[meanVE_abs_spatialD, sdLocR_abs_spatialD, seLocR_abs_spatialD] = ...
    deal(NaN(lenC, lenP, lenM, ceil(lenD/2))); 
sign_cor = [-1,-1,-1,1,1,1,1;1,1,1,1,-1,-1,-1];
%1st: A; 2nd: V; positive means shifting to the direction of the other modality

for i = 1:lenM
    for j = 1:lenC
        for k = 1:lenP
            %initialize cells that store individual localization responses
            %for each audiovisual spatial discrepancy
            %locR = cell(1,lenD);
            for l = 1:length(V_loc)
                for m = 1:length(A_loc)
                    locR_ijk = squeeze(locRespMat(j,k,m,l,i,:))';
                    %index:   m = 1, m = 2, m = 3, m = 4 (check the previous section)
                    %  l = 1    4,     3,     2,     1
                    %  l = 2    5,     4,     3,     2
                    %  l = 3    6,     5,     4,     3
                    %  l = 1    7,     6,     5,     4
                    locR{j,k,i,ceil(lenD/2)+l-m} = [locR{j,k,i,ceil(lenD/2)+l-m}, ...
                        locR_ijk - meanLocR_uni_rep{i}(l,m)];
                end
            end
            
            %absolute spatial discrepancy of 0
            meanVE(j,k,i,:) = arrayfun(@(idx) mean(locR{j,k,i,idx}), 1:lenD);
            sdLocR_pooled(j,k,i,:) = arrayfun(@(idx) std(locR{j,k,i,idx}), 1:lenD);
            seLocR_pooled(j,k,i,:) = arrayfun(@(idx) std(locR{j,k,i,idx})/sqrt(length(locR{j,k,i,idx})), 1:lenD);
            %sum across spatial discrepancy that have the same absolute values
            meanVE_abs_spatialD(j,k,i,1) = mean(locR{j,k,i,ceil(lenD/2)}); 
            sdLocR_abs_spatialD(j,k,i,1) = std(locR{j,k,i,ceil(lenD/2)});
            seLocR_abs_spatialD(j,k,i,1) = sdLocR_abs_spatialD(j,k,i,1)/...
                                            sqrt(length(locR{j,k,i,ceil(lenD/2)}));            
            
            %absolute spatial discrepancy of 8-24
            for n = 2:ceil(lenD/2)
                lb_idx_AVd = ceil(lenD/2)-n+1; ub_idx_AVd = ceil(lenD/2)+n-1;
                VE_corrected = [locR{j,k,i,lb_idx_AVd}.*sign_cor(i,lb_idx_AVd),...
                    locR{j,k,i,ub_idx_AVd}.*sign_cor(i,ub_idx_AVd)];
                meanVE_abs_spatialD(j,k,i,n) = mean(VE_corrected);
                sdLocR_abs_spatialD(j,k,i,n) = std(VE_corrected);
                seLocR_abs_spatialD(j,k,i,n) = sdLocR_abs_spatialD(j,k,i,n)/sqrt(length(VE_corrected));
            end
        end
    end
end

%% compute proportion of reporting C1 as a function of absolute spatial discrepancy
meanPC1_abs_spatialD = mean_pC1_lenD(:,:,ceil(lenD/2):end);
meanPC1_abs_spatialD(:,:,2:end) = meanPC1_abs_spatialD(:,:,2:end)./2 + ...
                                  mean_pC1_lenD(:,:,floor(lenD/2):-1:1)./2;
SD_pC1_abs_spatialD = NaN(size(meanPC1_abs_spatialD));
for i = 1:lenC
    for j = 1:lenP
        SD_pC1_abs_spatialD(i,j,:) = sqrt(squeeze(meanPC1_abs_spatialD(i,j,:)).*...
            (1-squeeze(meanPC1_abs_spatialD(i,j,:)))./numT_abs_AV');
    end
end
                    
%% plot the Ventriloquism effect as a function of spatial discrepancy
%info for plotting
y_bds_VE       = [AV_discrepancy(1)-2, AV_discrepancy(end)+2];
y_ticks_VE     = [AV_discrepancy(find(y_bds_VE(1) < AV_discrepancy,1):...
                             find(y_bds_VE(end) > AV_discrepancy,1,'last'))];
x_bds_VE       = [AV_discrepancy(1)-3, AV_discrepancy(end)+3];
x_ticks_VE     = [x_bds_VE(1), AV_discrepancy, x_bds_VE(end)];
cb             = [65,105,225]./255; cr = [0.85, 0.23, 0.28];
cMap_VE        = {min(cb.*2,1),cb; [247,191,190]./255, cr};
lgd_pos        = [0.3 0.8 0.05 0.1; 0.74 0.8 0.05 0.1];
y_lbl          = {'Auditory localization shifts (deg)',...
                  'Visual localization shifts (deg)'};

if bool_plt(3) == 1
    for i = 1:lenM
        figure
        for j = 1:lenC
            subplot(1,lenM,j)
            addBackground(x_bds_VE, y_bds_VE, x_ticks_VE, [y_bds_VE(1),...
                y_ticks_VE,y_bds_VE(end)])
            if i == 1
                plt3 = plot([AV_discrepancy(1),AV_discrepancy(end)],...
                    [AV_discrepancy(1),AV_discrepancy(end)],'k--',...
                    'Color',ones(1,3).*0.5,'lineWidth',lw); hold on;
            else
                plt3 = plot([AV_discrepancy(1),AV_discrepancy(end)],...
                    [AV_discrepancy(end),AV_discrepancy(1)],'k--',...
                    'Color',ones(1,3).*0.5,'lineWidth',lw); hold on;
            end
            plt4 = plot(AV_discrepancy, zeros(1,lenD), 'k:','Color',...
                ones(1,3).*0.5,'lineWidth',lw); hold on;
            for k = 1:lenP
                for n = 1:lenD
                    plt(k) = errorbar(AV_discrepancy(n)+jitter_prepost(k),...
                        meanVE(j,k,i,n), 2.*seLocR_pooled(j,k,i,n), '-o',...
                        'MarkerSize',sqrt(numT_AV(n)).*2, 'Color',...
                        cMap_VE{i,k},'MarkerFaceColor',cMap_VE{i,k},...
                        'MarkerEdgeColor',cMap_VE{i,k},'lineWidth',lw,...
                        'CapSize',15); hold on;
                end
                plot(AV_discrepancy+jitter_prepost(k), squeeze(meanVE(j,k,i,:)),...
                    '-', 'lineWidth',lw,'Color',cMap_VE{i,k}); hold on
            end
            text(x_bds_VE(1) + 7, y_bds_VE(1)+3, ['Condition: ', cond{j}],...
                'FontSize',fs_lgds); hold off;box off;
            %add legends
            if j == 1; text(x_bds_VE(1) + 0.2, y_bds_VE(end)-1, subjI,'FontSize',fs_lgds); end
            legend([plt(1) plt(2)], lgd_cond_phase,'Location','northeast',... %'Position',lgd_pos(j,:),...
                'FontSize',fs_lgds); legend boxoff;
            xticks(AV_discrepancy); xlim(x_bds_VE);
            xlabel('Spatial discrepancy (V-A, deg)');
            yticks(y_ticks_VE); ylim(y_bds_VE);
            ylabel(sprintf(y_lbl{i})); set(gca,'FontSize',fs_lbls);
        end
        %best size
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.75, 0.60]);
        set(gcf,'PaperUnits','centimeters','PaperSize',[50 25]);
%        saveas(gcf, ['VE_',modality{i}, 'shifts_', subjI], 'pdf');
        
%         %for manuscript
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.42, 0.40]);
%         set(gcf,'PaperUnits','centimeters','PaperSize',[40 20]);
%         saveas(gcf, ['VE_',modality{i}, 'shifts_', subjI,'_'], 'pdf');
    end
end

%% plot the Ventriloquism effect as a function of absolute spatial discrepancy
y_bds_VE_abs   = [-8, abs_AV_discrepancy(end)+2];
y_ticks_VE_abs = abs_AV_discrepancy;
x_bds_VE_abs   = [-3, abs_AV_discrepancy(end)+3];
x_ticks_VE_abs = [-3,abs_AV_discrepancy, abs_AV_discrepancy(end)+3];

if bool_plt(3) == 1
    for i = 1:lenM
        figure
        for j = 1:lenC
            subplot(1,lenM,j)
            addBackground(x_bds_VE_abs, y_bds_VE_abs, x_ticks_VE_abs, [y_bds_VE_abs(1),...
                y_ticks_VE_abs,y_bds_VE_abs(end)])
            if i == 1
                plt3 = plot(abs_AV_discrepancy,abs_AV_discrepancy,'k--',...
                    'Color',ones(1,3).*0.5,'lineWidth',lw); hold on;
            else
                plt3 = plot(abs_AV_discrepancy,abs_AV_discrepancy,'k--',...
                    'Color',ones(1,3).*0.5,'lineWidth',lw); hold on;
            end
            plt4 = plot(abs_AV_discrepancy, zeros(1,lenD_abs), 'k:','Color',...
                ones(1,3).*0.5,'lineWidth',lw); hold on;
            for k = 1:lenP
                for n = 1:lenD_abs
                    plt(k) = errorbar(abs_AV_discrepancy(n)+jitter_prepost(k),...
                        meanVE_abs_spatialD(j,k,i,n), 2.*seLocR_abs_spatialD(j,k,i,n),...
                        '-o','MarkerSize',sqrt(numT_abs_AV(n)).*2, 'Color',cMap_VE{i,k},...
                        'MarkerFaceColor',cMap_VE{i,k}, 'MarkerEdgeColor',...
                        cMap_VE{i,k},'lineWidth',lw,'CapSize',15); hold on;
                end
                plot(abs_AV_discrepancy+jitter_prepost(k), ...
                    squeeze(meanVE_abs_spatialD(j,k,i,:)),...
                    '-', 'lineWidth',lw,'Color',cMap_VE{i,k}); hold on
            end
            text(x_bds_VE_abs(1) + 7, y_bds_VE_abs(1)+3, ['Condition: ', cond{j}],...
                'FontSize',fs_lgds); hold off;box off;
            %add legends
            if j == 1; text(x_bds_VE_abs(1) + 0.2, y_bds_VE_abs(end)-1, subjI,'FontSize',fs_lgds); end
            legend([plt(1) plt(2)], lgd_cond_phase,'Location','northeast',... %'Position',lgd_pos(j,:),...
                'FontSize',fs_lgds); legend boxoff;
            xticks(abs_AV_discrepancy); xlim(x_bds_VE_abs);
            xlabel('Absolute spatial discrepancy (deg)');
            yticks(y_ticks_VE_abs); ylim(y_bds_VE_abs);
            ylabel(sprintf(y_lbl{i})); set(gca,'FontSize',fs_lbls);
        end
        %best size
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.65, 0.60]);
        set(gcf,'PaperUnits','centimeters','PaperSize',[50 25]); 
        saveas(gcf, ['VE_',modality{i}, 'shifts_abs_spatialD_', subjI], 'pdf');
    end
end

%% plot the unity judgment as a function of absolute spatial discrepancy
if bool_plt(1) == 1
    figure
    for i = 1:lenC
        subplot(1,2,i)
        addBackground(x_bds_VE_abs, y_bds, x_ticks_VE_abs, [y_bds(1), y_ticks, y_bds(end)])
        for j = 1:lenP
            for m = 1:lenD_abs
                plt(j) = errorbar(abs_AV_discrepancy(m)+jitter_prepost(j), ...
                    meanPC1_abs_spatialD(i,j,m), SD_pC1_abs_spatialD(i,j,m),'-o',...
                    'MarkerSize',sqrt(numT_abs_AV(m)).*2,'Color',cMap_unity(j,:),...
                    'MarkerFaceColor',cMap_unity(j,:),'MarkerEdgeColor',...
                    cMap_unity(j,:),'lineWidth',lw,'CapSize',15); hold on
            end
            plot(abs_AV_discrepancy+jitter_prepost(j), squeeze(meanPC1_abs_spatialD(i,j,:)),...
                '-','lineWidth',lw,'Color',cMap_unity(j,:)); hold on;
        end
        %add legends
        if i == 1;text(x_bds_VE_abs(1) + 0.2, 1.02, subjI,'FontSize',fs_lgds); hold on;end
        text(x_bds_VE_abs(1) + 5, 0.01, ['Condition: ', cond{i}], 'FontSize',fs_lgds); hold on
        legend([plt(1) plt(2)], lgd_cond_phase,'Location', 'northeast',...
            'FontSize',fs_lgds); legend boxoff;
        xticks(abs_AV_discrepancy); xlim(x_bds_VE_abs); xlabel('Absolute spatial discrepancy (deg)'); 
        ylabel(sprintf(['The probability \n of reporting a common cause'])); 
        yticks(y_ticks);ylim(y_bds);
        set(gca,'FontSize',fs_lbls);
    end
    %best size
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.65, 0.60]);
    set(gcf,'PaperUnits','centimeters','PaperSize',[50 25]);
    saveas(gcf, ['UnityJdg_btwPrePost_abs_spatialD_',subjI], 'pdf'); 
end

%% save data
BimodalLocalization_data = {subjN, unityMat, pC1_resp, locRespMat};
if bool_save(1) == 1
    save(['BimodalLocalization_sub', num2str(subjN), '_dataSummary.mat'],...
        'BimodalLocalization_data');
end

