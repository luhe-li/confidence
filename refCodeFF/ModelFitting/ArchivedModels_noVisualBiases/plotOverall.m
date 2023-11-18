%% load data
%best model (look-up table):
%'3-PW':full-MA-MAP,       '4-SW':full-MA-posterior,  '5-HL':null-MS-posterior,
%'6-YZ':full-MA-MAP,       '8-NH':samePC1pre-MA-MAP,  '9-ZZ':full-MA-MAP,
%'11-BB':full-MS-posterior,'12-ZY':full-MA-MAP,       '13-MR':null-MA-measurements,
%'18-ZL':samePC1pre-MA-MAP,'19-RE':full-MS-posterior, '20-MM':full-MA-MAP,
%'21-LL':samePC1pre-MA-measurements,                  '22-CS':null-MA-MAP,        
%'23-JH':samePC1pre-MA-MAP,'24-MD':samePC1pre-MA-MAP, '25-HHL':null-MA-MAP,
clear all; close all; clc

%subject number (outlier participants: 15, 16, 17)
subjN_dict = [   3,   4,   5,   6,   8,   9,  11,  12,  13,...
                18,  19,  20,  21,  22,  23,  24, 25];
%subject initial (outlier participants: 'AD', 'SM', 'SX'
subjI_dict = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR',...
              'ZL','RE','MM','LL','CS','JH','MD','HHL'};
%models that assume different plasticity of the common-cause prior
modelTypes_dict = {'samePC1pre', 'fullModel','null'};
%decision strategy
strategy_dict   = {'strategyMAP_MA_strategyUnity_posteriorC1',...
                   'strategyMAP_MA_strategyUnity_measurements',...
                   'strategyMAP_MA_strategyUnity_MAP',...
                   'strategyMAP_MS_strategyUnity_posteriorC1',...
                   'strategyMAP_MS_strategyUnity_measurements'};

%create popout window
prompt     = {'Subject ID (0: all subjects):',['Model of pCommon ',...
    '(1:samePC1pre, 2:full, 3:null, 0:all models)'], ['Model of strategy ',...
    '(1:MA-posteriorC1, 2:MA-measurements, 3:MA-MAP, 4: MS-posteriorC1,',...
    '5:MS-measurements, 0:all models)'], 'Fig.1: unity judgments (1:yes, 0:no)',...
    'Fig.2: localization responses (pre- & post-learning)',...
    'Fig.3: ventriloquism effects', 'Fig.4: delta AIC','Fig.5: #winnings',...
    'Fig.6: delta pCommon', ['Fig.7: Divide participants to',...
    'groups and plot separately'],'Save model comparison results'};
dlgtitle   = 'Input';
dims       = [1 100];
definput   = {'3','2','3','1','0','1','0','0','0','0','0'};
%to see group-level results
%definput  = {'0','0','0','0','0','0','1','1','1','1','0'}; 
answer     = inputdlg(prompt,dlgtitle,dims,definput);
bool_plot  = arrayfun(@(idx) str2double(answer(idx)), 4:10);  
bool_save  = zeros(1,7); %ones(1,7);
bool_saveD = str2double(answer(end));
%selected subjects   
if str2double(answer(1)) ~= 0 
    subjNs     = [str2double(answer(1))];
    subjIs     = {subjI_dict{subjNs==subjN_dict}};
else; subjNs = subjN_dict; subjIs = subjI_dict; 
end
    
%selected models
if str2double(answer(2)) ~= 0
    modelTypes = {modelTypes_dict{str2double(answer(2))}};
else; modelTypes = modelTypes_dict;
end

%selected strategies
if str2double(answer(3)) ~= 0 
    strategy = {strategy_dict{str2double(answer(3))}};
else; strategy = strategy_dict;
end
    
cond       = {'cong', 'incong'};
phase      = {'pre', 'post'};
modality   = {'Auditory','Visual'};
bool_data  = [1,1]; %including unity judgments, and localization responses
idx_mat    = {[8,9;8,10],[8,9;10,11],[8,8;8,8]}; %order is the same as modelTypes

%% Figs. 1-3: unity judgment, localization responses, ventriloquism effects
global sN sI bool_inc lenC lenP lenM lenMT lenDS lenS lenD
lenC  = length(cond); 
lenP  = length(phase); 
lenM  = length(modality);
lenMT = length(modelTypes);
lenDS = length(strategy);
lenS  = 4; %4 stimulus location (A or V)
lenD  = 7; %7 different spatial discrepancy
nS                 = length(subjNs);
[minNLL, deltaNLL] = deal(NaN(nS, lenMT, lenDS));
estimatedP         = cell(1,nS);
pC1_all            = NaN(nS, lenMT, lenDS, lenC, lenP);
for i = 1:nS
    sN = subjNs(i); sI = subjIs{i}; bool_inc = bool_data;
    [estimatedP{i}, minNLL(i,:,:)] = plotModelFits(modelTypes, strategy,...
        bool_plot(1:3), bool_save(1:3));
    for j = 1:lenMT
        for k = 1:lenDS
            pC1_all(i,j,k,1,:) = estimatedP{i}{j,k}(idx_mat{j}(1,:)); %cong condition
            pC1_all(i,j,k,2,:) = estimatedP{i}{j,k}(idx_mat{j}(2,:)); %incong condition
        end
    end
    deltaNLL(i,:,:) = minNLL(i,:,:) - min(min(minNLL(i,:,:)));
end

%% Fig.4: delta AIC
%We only need to run this section if all the 15 models are selected at the
%beginning
if str2double(answer(2)) == 0  && str2double(answer(3)) == 0
    green           = [76,153,0]./255;
    colorMap        = [linspace(green(1),1,255)', linspace(green(2),1,255)',...
                        linspace(green(3),1,255)'];
    getAIC          = @(nLL, k) 2.*nLL+ 2.*k;
    %getDeviance     = @(nLL_simple, nLL_complex) 2.*(nLL_simple - nLL_complex);

    %numParams is the total number of free parameters (note sigma_r is not
    %counted, and it's fine because the differences remain the same)
    numParams       = [12,13,13,12,13; 13,14,14,13,14; 10,11,11,10,11]; 
    %numTotalTrials  = 320 + 240 + 320*2*2; 
    %initialize
    [AIC, deltaAIC] = deal(NaN(nS, lenMT, lenDS));
    subjI_bestM     = cell(1,nS);
    subjI_bestM_sub = NaN(nS, 2); %subscripts

    %x and y tick labels
    xTicksL         = {'MA-Posterior-based','MA-Measurement-based',...
                       'MA-Estimate-based', 'MS-Posterior-based',...
                       'MS-Measurement-based'};
    xTicksL_short   = {'MA-PB','MA-MB', 'MA-EB', 'MS-PB','MS-MB'};
    yTicksL         = {'High-plasticity-short-lasting','High-plasticity-long-lasting', ...
                       'No-plasticity'};
    yTicksL_short   = {'HP-SL','HP-LL','NP'};

    for i = 1:nS
        %calculate delta AIC for each subject
        AIC(i,:,:)           = getAIC(squeeze(minNLL(i,:,:)), numParams);
        deltaAIC(i,:,:)      = AIC(i,:,:) - min(min(AIC(i,:,:)));
        %find the best model (the model that corresponds to 0 delta AIC)
        idx_bestM            = find(squeeze(deltaAIC(i,:,:)) == 0);
        [row, col]           = ind2sub([lenMT, lenDS], idx_bestM);
        subjI_bestM{i}       = [subjIs{i},'-', modelTypes{row},'-', strategy{col}];
        subjI_bestM_sub(i,:) = [row, col];
        disp([subjI_bestM{i}, ', nLL (AIC = 0): ', num2str(minNLL(i,row,col))]); %display it

        if bool_plot(4)
            figure
            heatmap(round(squeeze(deltaAIC(i,:,:)),1),'ColorbarVisible', 'on',...
                'YLabel', sprintf('Assumption about the\ncommon-cause prior'),...
                'XLabel', sprintf('Assumption about the decision strategy'), 'YData',...
                 yTicksL_short,'XData', xTicksL_short, 'GridVisible', 'off',...
                 'Colormap',colorMap); 
            caxis([0, 10]); set(gca,'FontSize',22); 
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.85, 0.45]);
            set(gcf,'PaperUnits','centimeters','PaperSize',[60 20]);
            if bool_save(4); saveas(gcf, ['ModelComparison_deltaAIC_',subjIs{i}], 'pdf'); end
        end
    end
end

%% Fig.5: #times a model wins
%We only need to run this section if all the participants and all the 
%models are selected at the beginning
if str2double(answer(1)) == 0 && str2double(answer(2)) == 0  && ...
        str2double(answer(3)) == 0 && bool_plot(5)
    %liberal criterion: 2; conservative criterion: 1e-4
    AIC_thres   = 1e-4; 
    counts_M    = squeeze(sum(deltaAIC < AIC_thres,1));
    org         = [255,140,0]./255;
    colorMap    = [linspace(1,org(1),255)',linspace(1,org(2),255)',...
                   linspace(1,org(3),255)'];

    figure
    heatmap(counts_M,'ColorbarVisible', 'on','YLabel', ...
        sprintf('Assumption about the\ncommon-cause prior'),'XLabel', ...
        sprintf('Assumption about the decision strategy'), 'YData',...
        yTicksL_short, 'XData', xTicksL_short,'GridVisible','off',...
        'Colormap',colorMap); 
    caxis([0, nS]); set(gca,'FontSize',22); 
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.85, 0.45]);
    set(gcf,'PaperUnits','centimeters','PaperSize',[60 20]);
    if bool_save(5)
        saveas(gcf, ['ModelComparison_counts_bestModel_DeltaAICthres',...
            num2str(AIC_thres)], 'pdf'); 
    end
    
    %count the number of times the full, restricted and the null model wins
    winCounts1 = sum(counts_M,2);
    winCounts2 = sum(counts_M,1);

    figure
    bar(1:3,winCounts1, 'FaceColor', [0.8,0.8,0.8]); ylim([0,20]); hold on
    for i = 1:lenMT
        text(i-1+0.95, winCounts1(i) + 0.8, num2str(winCounts1(i)), 'fontSize', 25);
    end
    hold off; box off; ylabel('Times of winning'); xticks(1:3); 
    xticklabels(yTicksL_short); yticks([]); xlim([0.25, 3.75]); 
    title(['$\Delta AIC = $', num2str(round(AIC_thres))],'Interpreter','latex'); 
    set(gca,'FontSize',22);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.45, 0.35]);
    set(gcf,'PaperUnits','centimeters','PaperSize',[30 14]);
    if bool_save(5)
        saveas(gcf, ['winCounts_ModelPC1_AICthres_', num2str(round(AIC_thres))], 'pdf');
    end

    figure
    bar(1:lenDS,winCounts2, 'FaceColor', [0.8,0.8,0.8]); ylim([0,20]); hold on
    for i = 1:lenDS
        text(i-1+0.95, winCounts2(i) + 1, num2str(winCounts2(i)), 'fontSize', 25);
    end
    hold off; box off; ylabel('Times of winning'); xticks(1:lenDS); 
    xticklabels(xTicksL_short); yticks([]); xlim([0.25, 5.75]); ylim([0,20]);
    title(['$\Delta AIC = $', num2str(round(AIC_thres))],'Interpreter','latex'); 
    set(gca,'FontSize',22);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.35]);
    set(gcf,'PaperUnits','centimeters','PaperSize',[60 14]);
    if bool_save(5)
        saveas(gcf, ['winCounts_Strategy_AICthres_', num2str(round(AIC_thres))], 'pdf');
    end
end

%% Fig.6: plot pC1 (fitted by the best model)
if str2double(answer(1)) == 0 && bool_plot(6)
    nBtst  = 100; %we bootstrap 100 times
    CI_bds = [0.025, 0.975];
    %first initialize matrices and cells
    %pC1_btst: each cell stores nBtst x #param estimated parameters
    %          (based on what the best-fitting parameter is)
    %pC1_btst_lb,ub: model estimates that correspond to 2.5% and 97.5% percentiles
    %diff_pC1_btst: the difference between pCommon_post and pCommon_pre
    %diff_pC1_btst_lb,ub: the diff values that correspond to 2.5% and 97.5% percentiles
    pC1_btst                             = cell(1,nS);
    [pC1_btst_lb,pC1_btst_ub]            = deal(NaN(nS, lenC, lenP));
    diff_pC1_btst                        = NaN(nS,nBtst,lenC);
    [diff_pC1_btst_lb, diff_pC1_btst_ub] = deal(NaN(nS, lenC));
    
    nLL_btst = NaN(nS, nBtst);
    [nLL_btst_lb, nLL_btst_ub] = deal(NaN(1,nS));
    %for each participant
    for i = 1:nS
        disp(i)
        %again load the best-fitting model
        bestM         = subjI_bestM{i};
        seg_idx       = find(bestM == '-'); 
        Model_name    = bestM((seg_idx(1)+1):((seg_idx(end)-1)));
        strategy_name = bestM((seg_idx(end)+1):end);
        %load the bootstrap file
        C_btst        = load(['ModelFitting_updatePrior_', Model_name,'_',...
                        strategy_name, '_btst_sub', num2str(subjNs(i)),'.mat'],...
                        'ModelFitting');
        %initialize
        pC1_btst{i}   = NaN(nBtst, length(C_btst.ModelFitting{end}{1}{end}));
        %get bootstrapped model estimates
        for j = 1:nBtst; pC1_btst{i}(j,:) = C_btst.ModelFitting{end}{j}{end}; end
        %the indices that correspond the estimated pCommon are different
        %for different models
        switch Model_name
            case 'samePC1pre';  idx = idx_mat{1};
            case 'fullModel';   idx = idx_mat{2};
            case 'null';        idx = idx_mat{3};
        end
        %initialize a cell that stores the best-estimated pC1 separately
        %for the pre- and the post-learning phase
        pC1_pre_post_btst = cell(1,lenP);
        for j = 1:lenP
            pC1_pre_post_btst{j} = [sort(pC1_btst{i}(:,idx(1,j))),...
                                    sort(pC1_btst{i}(:,idx(2,j)))];
            pC1_btst_lb(i,:,j)   = pC1_pre_post_btst{j}(ceil(nBtst*CI_bds(1)),:);
            pC1_btst_ub(i,:,j)   = pC1_pre_post_btst{j}(ceil(nBtst*CI_bds(2)),:);
        end
        %calculate the difference of best-estimated pC1 between the pre-
        %and the post-learning phases
        for k = 1:lenC
        	diff_pC1_btst(i,:,k) = sort(pC1_btst{i}(:,idx(k,2)) - pC1_btst{i}(:,idx(k,1)));
        end
        diff_pC1_btst_lb(i,:) = squeeze(diff_pC1_btst(i,ceil(nBtst*CI_bds(1)),:));
        diff_pC1_btst_ub(i,:) = squeeze(diff_pC1_btst(i,ceil(nBtst*CI_bds(2)),:));
        
        %get the nLL
        for j = 1:nBtst; nLL_btst(i,j) = C_btst.ModelFitting{end}{j}{3}; end
        nLL_btst_sorted = sort(nLL_btst(i,:));
        nLL_btst_lb(i) = nLL_btst_sorted(ceil(nBtst*CI_bds(1)));
        nLL_btst_ub(i) = nLL_btst_sorted(ceil(nBtst*CI_bds(2)));
    end
    
    pC1_bestM      = NaN(nS, lenC, lenP);
    diff_pC1_bestM = NaN(nS, lenC);
    choice_cmap    = 'uni'; %'uni' or 'rainbow'
    for i = 1:nS
        pC1_bestM(i,:,:) = squeeze(pC1_all(i, subjI_bestM_sub(i,1),...
            subjI_bestM_sub(i,2),:,:));
        diff_pC1_bestM(i,:) = pC1_bestM(i,:,2) - pC1_bestM(i,:,1);
    end
    plot_bestfitPC1(pC1_bestM,pC1_btst_lb,pC1_btst_ub,diff_pC1_bestM,...
        diff_pC1_btst_lb, diff_pC1_btst_ub, choice_cmap, subjIs, bool_save(6))
end

%% Fig. 7: group participants based on how the common-cause prior changes
%divide all the participants into three groups based on the patterns
if str2double(answer(1)) == 0 && bool_plot(7)
    grp  = {{'NH','LL','ZZ','MM'},{'ZL','JH','MD','YZ','BB','ZY','PW'},{'SW','RE'}};
    choice_cmap = 'rainbow';

    for i = 1:length(grp)
        %find the corresponding indices for those selected participants
        indices            = arrayfun(@(idx) find(strcmp(grp{i}{idx}, subjIs) == 1,1),...
                                1:length(grp{i}));
        pC1_bestM_i        = pC1_bestM(indices,:,:);
        diff_pC1_bestM_i   = diff_pC1_bestM(indices, :);
        pC1_btst_lb_i      = pC1_btst_lb(indices,:,:);
        pC1_btst_ub_i      = pC1_btst_ub(indices,:,:);
        diff_pC1_btst_lb_i = diff_pC1_btst_lb(indices,:,:);
        diff_pC1_btst_ub_i = diff_pC1_btst_ub(indices,:,:);
        plot_bestfitPC1(pC1_bestM_i, pC1_btst_lb_i, pC1_btst_ub_i, diff_pC1_bestM_i,...
            diff_pC1_btst_lb_i, diff_pC1_btst_ub_i, choice_cmap, grp{i}, bool_save(7), indices)
    end
end

%% save data
if str2double(answer(1)) == 0 && bool_saveD
    ModelComparison = {subjNs, subjIs, minNLL, deltaNLL, estimatedP, pC1_all,...
        AIC, deltaAIC,subjI_bestM};
    save('Results_modelComparison_overall.mat','ModelComparison');
end


