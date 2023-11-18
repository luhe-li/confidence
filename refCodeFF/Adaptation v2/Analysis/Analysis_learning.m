% load data from the learning phase
clear all; close all; clc
%subject number (outlier participants: 15, 16, 17)
subjN_dict = [   3,   4,   5,   6,   8,   9,  11,  12,  13, 18,  19,  20,...
                21,  22,  23,  24, 25];
lenS       = length(subjN_dict); 
lenM       = 2; %2 modalities: A, V
lenC       = 2; %2 conditions: congruent, incongruent
lenP       = 2; %2 phases: pre-learning, post-learning
lenT       = 2; %2 temporal offsets in the incongruent learning phase (-250, 250)
%subject initial (outlier participants: 'AD', 'SM', 'SX'
subjI_dict = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','ZL','RE','MM',...
              'LL','CS','JH','MD','HHL'};
%the order of condition was counterbalanced across participants
%1: congruent; 2: incongruent
order_dict = [ 2,1; 1,2; 2,1; 1,2; 2,1; 1,2; 2,1; 2,1; 1,2; 2,1; 1,2; 2,1;...
               2,1; 1,2; 2,1; 1,2; 2,1]; 

%% initialization
[locE_A, locE_A_1st, locE_A_2nd, locE_V, locE_V_1st, locE_V_2nd, ...
    STD_A, STD_A_1st, STD_A_2nd, STD_V, STD_V_1st, STD_V_2nd, ...
    prct_C1] = deal(NaN(lenS,lenC));
[STD_A_btst_CI, STD_A_btst_CI_1st, STD_A_btst_CI_2nd, STD_A_1,...
    STD_A_2, STD_A_3, STD_A_4, STD_V_1, STD_V_2, STD_V_3, STD_V_4,...
    STD_V_btst_CI, STD_V_btst_CI_1st, STD_V_btst_CI_2nd] = deal(NaN(lenS, lenC, 2));
[locE_A_incong2, locE_V_incong2, STD_A_incong2, STD_V_incong2] = deal(NaN(lenS, lenT)); 
[STD_A_incong2_btst_CI, STD_V_incong2_btst_CI] = deal(NaN(lenS, lenT, 2));

for ss = 1:lenS
    subjN     = subjN_dict(ss);
    subjI     = subjI_dict{ss};
    order_ses = order_dict(ss,:);
    % load data from unimodal spatial localization task
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'Unimodal localization v3/Data/', subjI]));
    E = load(['UnimodalLocalization_dataSummary_S', num2str(subjN),'.mat'],...
        'data_summary');
    mean_locResps_A = E.data_summary{2}{1};
    mean_locResps_V = E.data_summary{2}{2};

    % load files (congruent condition)
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'Adaptation v2/Data/', subjI]));
    C             = load(['Adaptation_congruent_sub', num2str(subjN),'_session',...
                        num2str(order_ses(1)),'.mat'], 'Adaptation_data');
    %1st row: AV loc (-12, -4, 4, 12), 2nd row: localization modality (1:A; 2:V)
    AVpairs_order{1} = C.Adaptation_data{end-4}.AVpairs_order;
    resp_modality{1} = C.Adaptation_data{end-4}.localize_modality;
    bool_unityJdg{1} = C.Adaptation_data{end-4}.bool_unityReport; %1:YES
    sLocs_V{1}       = C.Adaptation_data{3}.arrangedLocs_deg;
    sLocs_A{1}       = C.Adaptation_data{end-1}.arrangedLocs_deg;
    locResps{1}      = C.Adaptation_data{end}.localization;
    unityJdgs{1}     = C.Adaptation_data{end}.unity;
    prct_C1(ss,1)    = sum(unityJdgs{1} == 1)/sum(~isnan(unityJdgs{1}));

    %load files (incongruent condition)
    D = load(['Adaptation_incongruent_sub', num2str(subjN), '_session',...
        num2str(order_ses(2)), '.mat'], 'Adaptation_incongruent_data');
    AVpairs_order{2} = D.Adaptation_incongruent_data{end-4}.AVpairs_order;
    timing_V_incong  = D.Adaptation_incongruent_data{3}.timing_relative; 
    resp_modality{2} = D.Adaptation_incongruent_data{end-4}.localize_modality;
    bool_unityJdg{2} = D.Adaptation_incongruent_data{end-4}.bool_unityReport; %1:YES
    sLocs_V{2}       = D.Adaptation_incongruent_data{3}.arrangedLocs_deg;
    sLocs_A{2}       = D.Adaptation_incongruent_data{end-1}.arrangedLocs_deg;
    locResps{2}      = D.Adaptation_incongruent_data{end}.localization;
    unityJdgs{2}     = D.Adaptation_incongruent_data{end}.unity;
    prct_C1(ss,2)    = sum(unityJdgs{2} == 1)/sum(~isnan(unityJdgs{2}));

    %compute localization errors
    for i = 1:lenC
        %all data are selected
        [locE_A(ss,i), STD_A(ss,i), ~, STD_A_btst_CI(ss,i,:), ...
            locE_V(ss,i), STD_V(ss,i), ~, STD_V_btst_CI(ss,i,:)] = ...
            computeLocErrors(sLocs_A{i}, sLocs_V{i}, locResps{i}, ...
            resp_modality{i}, mean_locResps_A, mean_locResps_V);
        
        %we can also split the data into two halves, and see if sigma is
        %significantly different across the two halves
        [locE_A_1st(ss,i), STD_A_1st(ss,i), ~, STD_A_btst_CI_1st(ss,i,:), ...
            locE_V_1st(ss,i), STD_V_1st(ss,i), ~, STD_V_btst_CI_1st(ss,i,:)] = ...
            computeLocErrors(sLocs_A{i}, sLocs_V{i}, locResps{i}, ...
            resp_modality{i}, mean_locResps_A, mean_locResps_V, [],'1stHalf');
        
        [locE_A_2nd(ss,i), STD_A_2nd(ss,i), ~, STD_A_btst_CI_2nd(ss,i,:), ...
            locE_V_2nd(ss,i), STD_V_2nd(ss,i), ~, STD_V_btst_CI_2nd(ss,i,:)] = ...
            computeLocErrors(sLocs_A{i}, sLocs_V{i}, locResps{i}, ...
            resp_modality{i}, mean_locResps_A, mean_locResps_V, [],'2ndHalf');
        
        %we can also split the data into four quaters, and see if sigma is
        %significantly different across blocks
        [~, STD_A_1(ss,i), ~, ~, ~, STD_V_1(ss,i), ~, ~] = ...
            computeLocErrors(sLocs_A{i}, sLocs_V{i}, locResps{i}, ...
            resp_modality{i}, mean_locResps_A, mean_locResps_V, [],'1stQuarter');
        [~, STD_A_2(ss,i), ~, ~, ~, STD_V_2(ss,i), ~, ~] = ...
            computeLocErrors(sLocs_A{i}, sLocs_V{i}, locResps{i}, ...
            resp_modality{i}, mean_locResps_A, mean_locResps_V, [],'2ndQuarter');
        [~, STD_A_3(ss,i), ~, ~, ~, STD_V_3(ss,i), ~, ~] = ...
            computeLocErrors(sLocs_A{i}, sLocs_V{i}, locResps{i}, ...
            resp_modality{i}, mean_locResps_A, mean_locResps_V, [],'3rdQuarter');
        [~, STD_A_4(ss,i), ~, ~, ~, STD_V_4(ss,i), ~, ~] = ...
            computeLocErrors(sLocs_A{i}, sLocs_V{i}, locResps{i}, ...
            resp_modality{i}, mean_locResps_A, mean_locResps_V, [],'4thQuarter');
    end
    
    %split data based on the relative timing of V in the incongruent
    %condition
    timingOffsets = unique(timing_V_incong);
    idx_V_order   = {find(timing_V_incong < 0), find(timing_V_incong > 0)};
    for o = 1:length(timingOffsets)
        sLocs_A_slc       = sLocs_A{2}(idx_V_order{o});
        sLocs_V_slc       = sLocs_V{2}(idx_V_order{o});
        locResps_slc      = locResps{2}(idx_V_order{o});
        resp_modality_slc = resp_modality{2}(idx_V_order{o});
        [locE_A_incong2(ss,o), STD_A_incong2(ss,o), ~, ...
            STD_A_incong2_btst_CI(ss,o,:), locE_V_incong2(ss,o),...
            STD_V_incong2(ss,o), ~, STD_V_incong2_btst_CI(ss,o,:)] = ...
            computeLocErrors(sLocs_A_slc, sLocs_V_slc, locResps_slc, ...
            resp_modality_slc, mean_locResps_A, mean_locResps_V);
    end
end

%% organize data
locE               = {locE_A, locE_V};
locE_2halves       = {{locE_A_1st, locE_V_1st},{locE_A_2nd, locE_V_2nd}};
locE_group         = {mean(locE_A), mean(locE_V)};
locE_2halves_group = {{mean(locE_A_1st), mean(locE_V_1st)},{mean(locE_A_2nd), mean(locE_V_2nd)}};
SEM_group          = {std(locE_A)./sqrt(lenS),std(locE_V)./sqrt(lenS)};
SEM_2halves_group  = {{std(locE_A_1st)./sqrt(lenS), std(locE_V_1st)./sqrt(lenS)},...
                      {std(locE_A_2nd)./sqrt(lenS), std(locE_V_2nd)./sqrt(lenS)}};

STD                   = {STD_A, STD_V};
STD_2halves           = {{STD_A_1st, STD_V_1st},{STD_A_2nd, STD_V_2nd}};
STD_group             = {mean(STD_A), mean(STD_V)};
STD_2halves_group     = {{mean(STD_A_1st), mean(STD_V_1st)},{mean(STD_A_2nd), mean(STD_V_2nd)}};
SEM_STD_group         = {std(STD_A)./sqrt(lenS), std(STD_V)./sqrt(lenS)};
SEM_2halves_STD_group = {{std(STD_A_1st)./sqrt(lenS), std(STD_V_1st)./sqrt(lenS)},...
                         {std(STD_A_2nd)./sqrt(lenS), std(STD_V_2nd)./sqrt(lenS)}}; 
STD_btst_CI           = {STD_A_btst_CI, STD_V_btst_CI};
STD_2halves_btst_CI   = {{STD_A_btst_CI_1st, STD_V_btst_CI_1st},...
                         {STD_A_btst_CI_2nd,STD_V_btst_CI_2nd}};

%for incongruent condition, we need to repeat this for the two split
%datasets
STD_incong2           = {STD_A_incong2, STD_V_incong2};
STD_incong2_group     = {mean(STD_A_incong2), mean(STD_V_incong2)};
SEM_STD_incong2_group = {std(STD_A_incong2)./sqrt(lenS), std(STD_V_incong2)./sqrt(lenS)};
STD_incong2_btst_CI   = {STD_A_incong2_btst_CI, STD_V_incong2_btst_CI};

%
STD_4quarters_group = {{mean(STD_A_1), mean(STD_V_1)},{mean(STD_A_2), mean(STD_V_2)},...
    {mean(STD_A_3), mean(STD_V_3)},{mean(STD_A_4), mean(STD_V_4)}};
SEM_4quarters_STD_group = {{std(STD_A_1)./sqrt(lenS), std(STD_V_1)./sqrt(lenS)},...
    {std(STD_A_2)./sqrt(lenS), std(STD_V_2)./sqrt(lenS)},...
    {std(STD_A_3)./sqrt(lenS), std(STD_V_3)./sqrt(lenS)},...
    {std(STD_A_4)./sqrt(lenS), std(STD_V_4)./sqrt(lenS)}};

%% plot the mean localization errors
cmap  = [88,137,225; 210,111,115]./255;
x_bds = {[-20,20], [-15,15]}; x_ticks = {-20:10:20, -15:5:15};
modalityName = {'auditory','visual'};

figure
for m = 1:lenM
    subplot(1,lenM,m)
    addBackground(x_bds{m}, x_bds{m}, x_ticks{m}, x_ticks{m});
    %identity line
    plot(x_bds{m}, x_bds{m}, 'k--','lineWidth',1.25); hold on
    %individual data
    errorbar(locE{m}(:,1), locE{m}(:,2), STD{m}(:,2), STD{m}(:,2),...
        STD{m}(:,1), STD{m}(:,1),'.','Color',cmap(m,:),'lineWidth',0.5); hold on
    h1 = scatter(locE{m}(:,1), locE{m}(:,2),100,'o','filled','MarkerFaceColor',...
        cmap(m,:), 'MarkerFaceAlpha',0.2, 'MarkerEdgeColor',cmap(m,:)); 
    %group data
    errorbar(locE_group{m}(1), locE_group{m}(2), SEM_group{m}(2), ...
        SEM_group{m}(2), SEM_group{m}(1), SEM_group{m}(1), '.', 'Color',...
        max(cmap(m,:).*0.5, 0),'lineWidth',2);
    h2 = scatter(locE_group{m}(1), locE_group{m}(2),200,'d','filled',...
        'MarkerFaceColor', max(cmap(m,:).*0.5, 0),'MarkerEdgeColor',...
        max(cmap(m,:).*0.5, 0));
    %identity line
    hold off; box off; xlim(x_bds{m});ylim(x_bds{m}); axis square; 
    xticks(x_ticks{m}); yticks(x_ticks{m})
    xlabel(sprintf(['Mean ',modalityName{m},' localization errors ',...
        '\n(congruent learning phase)'])); 
    ylabel(sprintf(['Mean ',modalityName{m},' localization errors ',...
        '\n(incongruent learning phase)']));
    legend([h1,h2], {'Individual data','Group data (N=17)'},'Location',...
        'southeast', 'FontSize',18); legend boxoff
    set(gca,'FontSize',20);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.65]);
set(gcf,'PaperUnits','centimeters','PaperSize',[40 20]);
%saveas(gcf, 'MeanLocErrors_cong_vs_incong_learning', 'pdf');

%% plot the STD (localization variability with 95% bootstrap CI)
cmap    = [88,137,225; 210,111,115]./255;
x_bds   = {[0, 16],[0, 12]};
x_ticks = {0:4:16, 0:3:12};
modalityName = {'auditory','visual'};

figure
for m = 1:lenM
    subplot(1,lenM,m)
    addBackground(x_bds{m}, x_bds{m}, x_ticks{m}, x_ticks{m});
    %individual data (1st column: congruent; 2nd column: incongruent)
    errorbar(STD{m}(:,1), STD{m}(:,2), STD{m}(:,2) - STD_btst_CI{m}(:,2,1),...
        STD_btst_CI{m}(:,2,2) - STD{m}(:,2), STD{m}(:,1) - STD_btst_CI{m}(:,1,1),...
        STD_btst_CI{m}(:,1,2) - STD{m}(:,1),'.','Color',cmap(m,:),...
        'lineWidth',1.25); hold on
    h1 = scatter(STD{m}(:,1), STD{m}(:,2),100,cmap(m,:),'filled',...
        'MarkerFaceColor', cmap(m,:), 'MarkerFaceAlpha',0.3, ...
        'MarkerEdgeColor', cmap(m,:)); 
    %group data
    errorbar(STD_group{m}(1), STD_group{m}(2), SEM_STD_group{m}(2),...
        SEM_STD_group{m}(2), SEM_STD_group{m}(1),SEM_STD_group{m}(1), '.',...
        'Color',max(cmap(m,:).*0.5,0), 'lineWidth',2); 
    h2 = scatter(STD_group{m}(1), STD_group{m}(2),200, max(cmap(m,:).*0.5,0),...
        'filled', 'Marker','d'); 
    %identity line
    plot(x_bds{m}, x_bds{m},'k--','lineWidth',1.25); xlim(x_bds{m}); ylim(x_bds{m});
    hold off; box off; xticks(x_ticks{m}); yticks(x_ticks{m});
    axis square; xlabel(sprintf(['STD of ',modalityName{m},...
        ' localization errors \n(congruent learning phase)'])); 
    ylabel(sprintf(['STD of ', modalityName{m}, ' localization errors ',...
        '\n(incongruent learning phase)']));
    legend([h1,h2], {'Individual data','Group data (N=17)'},'Location',...
        'southeast', 'FontSize',18); legend boxoff
    set(gca,'FontSize',20);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.65]);
set(gcf,'PaperUnits','centimeters','PaperSize',[40 20]);
%saveas(gcf, 'STD_cong_vs_incong_learning', 'pdf');

%% plot the STD (localization variability with 95% bootstrap CI) but here 
%we split the data into two halves, and see if the variability changes
%over time
markerStyle  = {'^','v'};

figure
for m = 1:lenM
    subplot(1,lenM,m)
    addBackground(x_bds{m}, x_bds{m}, x_ticks{m}, x_ticks{m});
    %individual data (1st column: congruent; 2nd column: incongruent)
    for f = 1:2
        errorbar(STD_2halves{f}{m}(:,1), STD_2halves{f}{m}(:,2),...
            STD_2halves{f}{m}(:,2) - STD_2halves_btst_CI{f}{m}(:,2,1),...
            STD_2halves_btst_CI{f}{m}(:,2,2) - STD_2halves{f}{m}(:,2),...
            STD_2halves{f}{m}(:,1) - STD_2halves_btst_CI{f}{m}(:,1,1),...
            STD_2halves_btst_CI{f}{m}(:,1,2) - STD_2halves{f}{m}(:,1),...
            '.','Color',cmap(m,:), 'lineWidth',0.75); hold on
        g(f) = scatter(STD_2halves{f}{m}(:,1), STD_2halves{f}{m}(:,2),100,...
            cmap(m,:),'filled','MarkerFaceColor', cmap(m,:), 'MarkerFaceAlpha',0.3, ...
            'MarkerEdgeColor', cmap(m,:),'Marker',markerStyle{f}); 
    end
    %group data
    for f = 1:2
        errorbar(STD_2halves_group{f}{m}(1), STD_2halves_group{f}{m}(2),...
            SEM_2halves_STD_group{f}{m}(2), SEM_2halves_STD_group{f}{m}(2),...
            SEM_2halves_STD_group{f}{m}(1),SEM_2halves_STD_group{f}{m}(1), '.',...
            'Color',max(cmap(m,:).*0.5,0), 'lineWidth',2); 
        scatter(STD_2halves_group{f}{m}(1), STD_2halves_group{f}{m}(2),200,...
            max(cmap(m,:).*0.5,0),'filled', 'Marker',markerStyle{f}); 
    end
    %identity line
    plot(x_bds{m}, x_bds{m},'k--','lineWidth',1.25); xlim(x_bds{m}); ylim(x_bds{m});
    hold off; box off; xticks(x_ticks{m}); yticks(x_ticks{m});
    axis square; xlabel(sprintf(['STD of ',modalityName{m},...
        ' localization errors \n(congruent learning phase)'])); 
    ylabel(sprintf(['STD of ', modalityName{m}, ' localization errors ',...
        '\n(incongruent learning phase)']));
    legend(g, {'1st half of the data', '2nd half of the data'}, 'Location',...
        'southeast','FontSize',18); legend boxoff
    set(gca,'FontSize',20);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.65]);
set(gcf,'PaperUnits','centimeters','PaperSize',[40 20]);
%saveas(gcf, 'STD_cong_vs_incong_learning_2halves', 'pdf');

%% plot the STD (localization variability with 95% bootstrap CI) but here 
%we split the data into four quarters, and see if the variability changes
%over time
markerStyle  = {'^','v','<','>'};

figure
for m = 1:lenM
    subplot(1,lenM,m)
    addBackground(x_bds{m}, x_bds{m}, x_ticks{m}, x_ticks{m});
    %group data
    for f = 1:4
        errorbar(STD_4quarters_group{f}{m}(1), STD_4quarters_group{f}{m}(2),...
            SEM_4quarters_STD_group{f}{m}(2), SEM_4quarters_STD_group{f}{m}(2),...
            SEM_4quarters_STD_group{f}{m}(1),SEM_4quarters_STD_group{f}{m}(1), '.',...
            'Color',max(cmap(m,:).*0.5,0), 'lineWidth',1); 
        g(f) = scatter(STD_4quarters_group{f}{m}(1), STD_4quarters_group{f}{m}(2),200,...
            max(cmap(m,:).*0.5,0),'filled', 'Marker',markerStyle{f},...
            'MarkerFaceAlpha',0.4); 
    end
    %identity line
    plot(x_bds{m}, x_bds{m},'k--','lineWidth',1.25); xlim(x_bds{m}); ylim(x_bds{m});
    hold off; box off; xticks(x_ticks{m}); yticks(x_ticks{m});
    axis square; xlabel(sprintf(['STD of ',modalityName{m},...
        ' localization errors \n(congruent learning phase)'])); 
    ylabel(sprintf(['STD of ', modalityName{m}, ' localization errors ',...
        '\n(incongruent learning phase)']));
    legend(g, {'1st quarter of the data', '2nd quarter of the data',...
        '3rd quarter of the data','4th quarter of the data'}, 'Location',...
        'southeast','FontSize',18); legend boxoff
    set(gca,'FontSize',20);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.65]);
set(gcf,'PaperUnits','centimeters','PaperSize',[40 20]);
%saveas(gcf, 'STD_cong_vs_incong_learning_4quarters', 'pdf');


%% plot the STD (localization variability with 95% bootstrap CI) but here 
%we split the data from the incongruent learning phase based on the
%temporal order between V and A
%auditory
cmap         = [88,137,225; 210,111,115]./255;
x_bds        = {[0, 16],[0, 12]};
x_ticks      = {0:4:16, 0:3:12};
modalityName = {'auditory','visual'};
markerStyle  = {'<','>'};

figure
for m = 1:lenM
    subplot(1,lenM,m)
    addBackground(x_bds{m}, x_bds{m}, x_ticks{m}, x_ticks{m});
    %individual data (1st column: congruent; 2nd column: incongruent)
    for i = 1:lenT
        errorbar(STD{m}(:,1), STD_incong2{m}(:,i), STD_incong2{m}(:,i) - ...
            STD_incong2_btst_CI{m}(:,i,1), STD_incong2_btst_CI{m}(:,i,2) - ...
            STD_incong2{m}(:,i), STD{m}(:,1) - STD_btst_CI{m}(:,1,1),...
            STD_btst_CI{m}(:,1,2) - STD{m}(:,1),'.','Color',cmap(m,:),...
            'lineWidth',0.75); hold on
        h(i) = scatter(STD{m}(:,1), STD_incong2{m}(:,i),100,cmap(m,:),'filled',...
            'MarkerFaceColor', cmap(m,:), 'MarkerFaceAlpha',0.3,...
            'MarkerEdgeColor', cmap(m,:),'Marker',markerStyle{i}); 
    end
    for i = 1:lenT
        %group data
        errorbar(STD_group{m}(1), STD_incong2_group{m}(i), ...
            SEM_STD_incong2_group{m}(i), SEM_STD_incong2_group{m}(i),...
            SEM_STD_group{m}(1),SEM_STD_group{m}(1), '.', 'Color',...
            max(cmap(m,:).*0.5,0), 'lineWidth',2); 
        scatter(STD_group{m}(1), STD_incong2_group{m}(i),200, ...
            max(cmap(m,:).*0.5,0), 'filled', 'Marker',markerStyle{i}); 
    end
    %identity line
    plot(x_bds{m}, x_bds{m},'k--','lineWidth',1.25); xlim(x_bds{m}); ylim(x_bds{m}); 
    hold off; box off; xticks(x_ticks{m}); yticks(x_ticks{m});
    axis square; xlabel(sprintf(['STD of ', modalityName{m},...
        ' localization errors \n(congruent learning phase)'])); 
    ylabel(sprintf(['STD of ',modalityName{m},' localization errors ',...
        '\n(incongruent learning phase)']));
    legend(h, {'$SOA(t_V - t_A) = -350 ms$', '$SOA(t_V - t_A) = 350 ms$'},...
        'interpreter','latex','Location','southeast','FontSize',18); legend boxoff
    set(gca,'FontSize',20);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.65]);
set(gcf,'PaperUnits','centimeters','PaperSize',[40 20]);
%saveas(gcf, 'STD_cong_vs_incong_learning_incong_split', 'pdf');

%% load estimated sigma_AV_A, sigma_AV_V for the pre- and the post-learning phases
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'ModelFitting/Fits/']));
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'Pointing task/Data/']));
G = load('Results_sigma_r.mat', 'sigma_r');
if lenS == length(G.sigma_r{3}); sigma_r = G.sigma_r{3}; end

addpath(genpath('/Users/hff/Desktop/NYU/Project2/Experiment code/ModelFitting/'));
F               = load('Results_modelComparison_overall.mat', 'ModelComparison');
estimatedP      = F.ModelComparison{5};
bestM           = F.ModelComparison{end};
pC1_all         = F.ModelComparison{end-3}; %size: lenS, lenMT, lenDS, lenC, lenP
%models that assume different plasticity of the common-cause prior
modelTypes_dict = {'samePC1pre', 'fullModel','null'};
%decision strategy
strategy_dict   = {'strategyMAP_MA_strategyUnity_posteriorC1',...
                   'strategyMAP_MA_strategyUnity_measurements',...
                   'strategyMAP_MA_strategyUnity_MAP',...
                   'strategyMAP_MS_strategyUnity_posteriorC1',...
                   'strategyMAP_MS_strategyUnity_measurements'};
[est_sigma_AV, STD_prepost] = deal(NaN(lenS, lenM)); %sigma_AV_A, sigma_AV_V
STD_btst_CI_prepost         = NaN(lenS, lenM, 2); %lower bound, upper bounds
pC1_bestM                   = NaN(lenS, lenC, lenP);
delta_pC1                   = NaN(lenS, lenC);
delta_ratio_pC1             = NaN(lenS, lenC);

for s = 1:lenS
    bestM_s           = bestM{s};
    seg_idx           = find(bestM_s == '-'); 
    Model_name        = bestM_s((seg_idx(1)+1):((seg_idx(end)-1)));
    strategy_name     = bestM_s((seg_idx(end)+1):end);
    Model_idx         = find(strncmp(Model_name, modelTypes_dict, 1)==1);
    strategy_idx      = find(strncmp(strategy_dict, strategy_name, 32)==1);
    est_sigma_AV(s,:) = estimatedP{s}{Model_idx, strategy_idx}(5:6);
    %when we take motor noise into account
    STD_prepost(s,:)  = sqrt(est_sigma_AV(s,:).^2 + sigma_r(s)^2);
    pC1_bestM(s,:,:)  = pC1_all(s, Model_idx, strategy_idx, :, :);
    delta_pC1(s,:)    = pC1_bestM(s,:,2) - pC1_bestM(s,:,1);
    delta_ratio_pC1(s,:)=delta_pC1(s,:)./pC1_bestM(s,:,1);
    
    %load bootstrap CI on estimated parameters
    G                 = load(['ModelFitting_updatePrior_', Model_name, '_',...
                        strategy_name,'_btst_sub',num2str(subjN_dict(s)),'.mat']);
    est_sigma_AV_btst = NaN(100,2);
    for b = 1:100; est_sigma_AV_btst(b,:) = G.ModelFitting{end}{b}{end}(5:6);end
    %sort them 
    est_sigma_AV_btst_sort     = sort(est_sigma_AV_btst);
    STD_btst_CI_prepost(s,:,1) = est_sigma_AV_btst_sort(3,:); %lower bounds
    STD_btst_CI_prepost(s,:,2) = est_sigma_AV_btst_sort(98,:); %upper bounds
end

%compute mean sigma_AV_r_prepost
STD_group_prepost = mean(STD_prepost);
SEM_group_prepost = std(STD_prepost)./sqrt(lenS);

%%
condName = {'congruent','incongruent'};
x_bds = [-3, 17]; x_ticks = 0:5:17;

figure
for p = 1:lenC
    subplot(1,lenC,p)
    addBackground(x_bds, x_bds, [x_bds(1), x_ticks, x_bds(end)],...
        [x_bds(1), x_ticks, x_bds(end)]);
    %identity line
    plot(x_bds, x_bds, 'k--','lineWidth',1.25); hold on
    for m = 1:lenM
        %individual data
        errorbar(STD_prepost(:,m), STD{m}(:,p), STD{m}(:,p) - ...
            STD_btst_CI{m}(:,p,1),STD_btst_CI{m}(:,p,2)-STD{m}(:,p), ...
            STD_prepost(:,m)-STD_btst_CI_prepost(:,m,1),...
            STD_btst_CI_prepost(:,m,2) - STD_prepost(:,m),'.',...
            'Color',cmap(m,:),'lineWidth',0.75); hold on
        h5(m) = scatter(STD_prepost(:,m), STD{m}(:,p), 100,cmap(m,:),...
            'filled', 'MarkerFaceAlpha',0.2,'MarkerEdgeColor', cmap(m,:),...
            'MarkerEdgeAlpha',1); hold on
        %group data
        errorbar(STD_group_prepost(m), STD_group{m}(p), SEM_group{m}(p),...
            SEM_group{m}(p), SEM_group_prepost(m), SEM_group_prepost(m),...
            '.','Color',max(cmap(m,:).*0.5,0),'lineWidth',2); 
        h6(m) = scatter(STD_group_prepost(m), STD_group{m}(p), 150,cmap(m,:),...
            'filled','Marker','d','MarkerFaceColor', max(cmap(m,:).*0.5,0),...
            'MarkerEdgeColor',max(cmap(m,:).*0.5,0));
    end
    hold off; box off; xlim(x_bds); ylim(x_bds); axis square;
    xlabel(sprintf('Localization noise \n(pre- & post-learning phases)'));
    ylabel(sprintf(['Localization noise \n(', condName{p},' learning phase)']));
    if p==1
        legend([h5(1),h6(1)], {'Individual data (auditory)','Group data (N=17)'},...
            'Location','southeast','FontSize',18); 
    else
        legend([h5(2),h6(2)], {'Individual data (visual)','Group data (N=17)'},...
            'Location','southeast','FontSize',18); 
    end
    legend boxoff; set(gca,'FontSize',20);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.65]);
set(gcf,'PaperUnits','centimeters','PaperSize',[40 20]);
%saveas(gcf, 'Sigma_preNpost_vs_learning', 'pdf');

%% save data
data_learning.locE                   = locE; 
data_learning.locVar                 = STD; 
data_learning.locVar_btst_CI         = STD_btst_CI; 
data_learning.locE_group             = locE_group;
data_learning.SEM_group              = SEM_group;

data_learning.locE_incong2           = locE_A_incong2; 
data_learning.STD_incong2            = STD_incong2;
data_learning.locVar_incong2_btst_CI = STD_incong2_btst_CI;
data_learning.STD_incong2_group      = STD_incong2_group;
data_learning.SEM_STD_incong2_group  = SEM_STD_incong2_group;

data = {subjN_dict, subjI_dict, data_learning, sigma_r};
save('Analysis_learningPhase.mat','data');

%% plot the changes in pC1 as a function of changes in sigma
% subjN_dict = [   3,   4,   5,   6,   8,   9,  11,  12,  13, 18,  19,  20,...
%                 21,  22,  23,  24, 25];
% subjI_dict = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','ZL','RE','MM',...
%               'LL','CS','JH','MD','HHL'};
% grp  = {{'NH','LL','ZZ','MM'},{'ZL','JH','MD','YZ','BB','ZY','PW'},{'SW','RE'}};


subjN_cong = {[   9,  20,  18,  23,   6,  12,   3],...
              [   8,  21,   9,  20,  24,   6,  11,  3,  4,  19]};
subjI_cong = {{'ZZ','MM','ZL','JH','YZ','ZY','PW'},...
              {'NH','LL','ZZ','MM','MD','YZ','BB','PW','SW','RE'}};

figure
for c = 1:lenC
    subplot(1,lenC,c)
    if c==1; addBackground([-1,0],[-0.4, 0.8],-1:0.5:0.5, -0.4:0.4:0.8);
    else; addBackground([0,3], [-0.3, 0.7], -1:1:3, -0.3:0.3:0.7); end
    
    idx_s = arrayfun(@(s) find(subjN_cong{c}(s) == subjN_dict), 1:length(subjN_cong{c}));
    for m = 1:lenM
        scatter((STD{m}(idx_s,c) - STD_prepost(idx_s,m))./STD_prepost(idx_s,m),...
            delta_pC1(idx_s,c), 150,cmap(m,:), 'filled','MarkerFaceColor',...
            cmap(m,:), 'MarkerFaceAlpha',0.4, 'MarkerEdgeColor', cmap(m,:)); hold on
    end
    plot([0,0],[-0.4, 0.8], 'k-'); plot([-1,3],[0,0],'k--'); hold off
    if c == 1; ylim([-0.4,0.8]); xlim([-1,0.5]);
    else; ylim([-0.3, 0.7]); xlim([-1,3]); end
    xlabel(sprintf('Localization noise of learning - of pre & post \n /localization noise of pre & post'));
    if c == 1; ylabel(sprintf('Changes in the common-cause prior\n after congruent learning phase'));
    else; ylabel(sprintf('Changes in the common-cause prior\n after incongruent learning phase')); end
    set(gca,'FontSize',20);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.8, 0.55]);
set(gcf,'PaperUnits','centimeters','PaperSize',[40 22]);
%saveas(gcf, 'deltaPC1_vs_deltaSigma', 'pdf');

%% plot changes of pCommon and the proportion of reporting C=1 during the learning phase
figure
for c = 1:lenC
    subplot(1,lenC,c)    
    idx_s = arrayfun(@(s) find(subjN_cong{c}(s) == subjN_dict), 1:length(subjN_cong{c}));
    scatter(prct_C1(idx_s,c), delta_pC1(idx_s,c), 150, ones(1,3).*0.5, ...
        'filled','MarkerFaceColor', ones(1,3).*0.5, 'MarkerFaceAlpha',0.4,...
        'MarkerEdgeColor', cmap(m,:)); hold on
    hold off
    if c == 1; ylim([-0.4,0.8]); xlim([0,1]);
    else; ylim([-0.3, 0.7]); xlim([0,1]); end
    xlabel(sprintf('Localization noise of learning - of pre & post \n /localization noise of pre & post'));
    if c == 1; ylabel(sprintf('Changes in the common-cause prior\n after congruent learning phase'));
    else; ylabel(sprintf('Changes in the common-cause prior\n after incongruent learning phase')); end
    set(gca,'FontSize',20);
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.8, 0.55]);
set(gcf,'PaperUnits','centimeters','PaperSize',[40 22]);





