%% subject info
clear all; close all; clc
subjNs   = [   3,   4,   5,   6,   8,   9,  11,  12,  13,  18,  19,  20,...
              21,  22,  23,  24,  25];  %15AD, 16SM, 17SX are identified as outliers
subjIs   = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','ZL','RE','MM',...
            'LL','CS','JH','MD','HHL'};

%% put all xlsx files together
varName1 = {'UnityJudgment','Condition','Phase','SpatialD','SpatialD_abs','SubjI'};
varName2 = {'LocShifts','LocShifts_corr','Condition','Phase','Modality',...
    'SpatialD','SpatialD_abs','SubjI'};
[sI_cell, cond_cell, phase_cell, modality_cell, spatialD_mat, phase_abs_AVd,...
    spatialD_abs_mat, unity_mat,locShifts_mat, locShifts_corr_mat, ...
    cond_abs_AVd, locResp_abs_AVd, modality_abs_AVd, sI_abs_AVd] = deal([]); %phase_abs_AVd, spatialD_abs
for i = 1:length(subjNs) 
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'Bimodal localization v2/Data/', subjIs{i}]));
    T1 = readtable(['unityJdg_binaryData_', subjIs{i}, '.xlsx']);
    unity_mat    = [unity_mat; T1.UnityJudgment]; 
    cond_cell    = [cond_cell; T1.Condition];
    phase_cell   = [phase_cell; T1.Phase];
    spatialD_mat = [spatialD_mat; T1.SpatialD];
    sI_cell      = [sI_cell; T1.SubjI];
    
    T2 = readtable(['locResp_demeanedData_', subjIs{i}, '.xlsx']);
    modality_cell = [modality_cell; T2.Modality];
    locShifts_mat   = [locShifts_mat; T2.LocResp];
    locShifts_corr_mat = [locShifts_corr_mat; T2.LocResp_crr];
    
    T3 = readtable(['locResp_demeanedData_abs_AVd_', subjIs{i}, '.xlsx']);
    modality_abs_AVd = [modality_abs_AVd; T3.Modality];
    cond_abs_AVd = [cond_abs_AVd; T3.Condition];
    phase_abs_AVd = [phase_abs_AVd; T3.Phase];
    %spatialD_abs = [spatialD_abs; T3.SpatialD];
    sI_abs_AVd = [sI_abs_AVd; T3.SubjI];
    locResp_abs_AVd =[locResp_abs_AVd; T3.LocResp_abs_spatialD];
end

%% overall unity judgment (raw data)
Table_unity = table(unity_mat, cond_cell, phase_cell, spatialD_mat,...
    abs(spatialD_mat), sI_cell, 'VariableNames', varName1);
fileName1 = 'unityJdg_binaryData_overall.xlsx';
writetable(Table_unity,fileName1,'Sheet','MyNewSheet','WriteVariableNames',true);
writetable(Table_unity,'unityJdg_binaryData_overall.txt');

%% overall demeaned localization responses
Table_locResp = table(locShifts_mat,locShifts_corr_mat, cond_cell, ...
    phase_cell, modality_cell, spatialD_mat, abs(spatialD_mat), sI_cell,...
    'VariableNames', varName2);
fileName2 = 'locShifts_indvdTrials_overall.xlsx';
writetable(Table_locResp,fileName2,'Sheet','MyNewSheet','WriteVariableNames',true);
writetable(Table_locResp,'locShifts_indvdTrials_overall.txt');

% %% average ventriloquism effects
% Table_locResp_abs_AVd = table(locResp_abs_AVd, cond_abs_AVd, phase_abs_AVd, ...
%     modality_abs_AVd,spatialD_abs, sI_abs_AVd,'VariableNames', varName2);
% fileName3 = 'locResp_demeanedData_abs_AVd_overall.xlsx';
% writetable(Table_locResp_abs_AVd,fileName3,'Sheet','MyNewSheet',...
%     'WriteVariableNames',true);

%% put all xlsx files together for percent of C=1 responses and VE
varName4 = {'prctC1','nTrials','nTrialsC1','Condition','Phase','SpatialD',...
    'SpatialD_abs','SubjI', 'subjN'};
varName5 = {'VE_Ashift','VE_Ashift_crr','nTrials','Condition','Phase','SpatialD',...
    'SpatialD_abs','SubjI', 'subjN'};
varName6 = {'VE_Vshift','VE_Vshift_crr','nTrials','Condition','Phase','SpatialD',...
    'SpatialD_abs','SubjI', 'subjN'};
[sI_cell,sN_mat,cond_cell, phase_cell, spatialD_mat, prctC1_mat, ...
    VE_Ashift_mat, nTrials,VE_Vshift_mat, VE_Ashift_crr, VE_Vshift_crr,...
    phase_abs_AVd, nTrials_mat] = deal([]);
lenC = 2; lenP = 2; lenD = 7; nT = repmat([1:1:4,3:-1:1]'.*20,[lenC*lenP,1]);
[crr_sign_A, crr_sign_V] = deal(ones(lenC*lenP*lenD, 1));
crr_sign_A(1:lenD:end) = -1; crr_sign_A(2:lenD:end) = -1; crr_sign_A(3:lenD:end) = -1;
crr_sign_V(5:lenD:end) = -1; crr_sign_V(6:lenD:end) = -1; crr_sign_V(7:lenD:end) = -1;

for i = 1:length(subjNs) 
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'Bimodal localization v2/Data/', subjIs{i}]));
    T4            = readtable(['prctC1_', subjIs{i}, '.xlsx']);
    prctC1_mat    = [prctC1_mat; T4.PercentC1]; 
    cond_cell     = [cond_cell; T4.Condition];
    nTrials_mat   = [nTrials_mat; nT];
    phase_cell    = [phase_cell; T4.Phase];
    spatialD_mat  = [spatialD_mat; T4.SpatialD];
    sN_mat        = [sN_mat; ones(length(T4.SubjI),1).*i];
    sI_cell       = [sI_cell; T4.SubjI];
    
    T5            = readtable(['VE_Ashift_', subjIs{i}, '.xlsx']);
    VE_Ashift_mat = [VE_Ashift_mat; T5.Ashift_VE];
    VE_Ashift_crr = [VE_Ashift_crr; T5.Ashift_VE.*crr_sign_A];
    
    T6            = readtable(['VE_Vshift_', subjIs{i}, '.xlsx']);
    VE_Vshift_mat = [VE_Vshift_mat; T6.Vshift_invVE];
    VE_Vshift_crr = [VE_Vshift_crr; T6.Vshift_invVE.*crr_sign_V];
end
spatialD_abs = abs(spatialD_mat);
nTrials_C1_mat = round(nTrials_mat.*prctC1_mat);

Table_prctC1 = table(prctC1_mat, nTrials_mat,nTrials_C1_mat, cond_cell,...
    phase_cell,spatialD_mat, spatialD_abs, sI_cell,sN_mat, 'VariableNames',...
    varName4);
writetable(Table_prctC1, 'prctC1_overall.xlsx','Sheet','MyNewSheet',...
    'WriteVariableNames',true);
writetable(Table_prctC1,'prctC1_overall.txt');

Table_Ashift = table(VE_Ashift_mat,VE_Ashift_crr, nTrials_mat./2,...
    cond_cell, phase_cell, spatialD_mat, spatialD_abs, sI_cell,sN_mat,...
    'VariableNames', varName5);
writetable(Table_Ashift,'Ashift_overall.xlsx','Sheet','MyNewSheet',...
    'WriteVariableNames',true);
writetable(Table_Ashift,'Ashift_overall.txt');

Table_Vshift = table(VE_Vshift_mat,VE_Vshift_crr, nTrials_mat./2,...
    cond_cell, phase_cell, spatialD_mat, spatialD_abs, sI_cell,sN_mat,...
    'VariableNames', varName6);
writetable(Table_Vshift,'Vshift_overall.xlsx','Sheet','MyNewSheet',...
    'WriteVariableNames',true);
writetable(Table_Vshift,'Vshift_overall.txt');

%% do t-tests on auditory and visual shifts
len_abs_AVd = 4;
locResp_abs_AVd_subj  = reshape(locResp_abs_AVd, ...
    [length(locResp_abs_AVd)/length(subjNs), length(subjNs)]);
idx_temp = reshape(1:(length(locResp_abs_AVd)/length(subjNs)),4,...
    (length(locResp_abs_AVd)/length(subjNs))/4);
idx_A = idx_temp(:,1:2:end);
idx_V = idx_temp(:,2:2:end);
locRespA_abs_AVd_subj = locResp_abs_AVd_subj(idx_A(:)',:);
locRespV_abs_AVd_subj = locResp_abs_AVd_subj(idx_V(:)',:);

%two-sample t test to examine whether auditory and visual ventriloquism
%effects are significantly different from each other
[~, P, STATS,d] = af_cohenD(mean(locRespA_abs_AVd_subj), mean(locRespV_abs_AVd_subj));
disp(STATS)
disp(P)
disp(d)


%one-sample t test to examine whether auditory ventriloquism effects are
%significantly different from 0
for i = 1:len_abs_AVd
    locRespA_abs_AVd_subj_i = mean(locRespA_abs_AVd_subj(i:len_abs_AVd:end,:),1);
    [~, P, STATS,d] = af_cohenD(locRespA_abs_AVd_subj_i);
    disp(['A, abs spatial discrepancy level: ', num2str(i)])
    disp(STATS)
    disp(P)
    disp(d)
end

%one-sample t test to examine whether visual ventriloquism effects are
%significantly different from 0
for i = 1:len_abs_AVd
    locRespV_abs_AVd_subj_i = mean(locRespV_abs_AVd_subj(i:len_abs_AVd:end,:),1);
    [~, P, STATS,d] = af_cohenD(locRespV_abs_AVd_subj_i);
    disp(['V, abs spatial discrepancy level: ', num2str(i)])
    disp(STATS)
    disp(P)
    disp(d)
end

%% do t-tests on spatial discrepancy
unique_absD = unique(spatialD_abs);
modality    = unique(modality_abs_AVd);
len_absD    = length(unique_absD);
lenM        = length(modality);
comp_idx    = [1,2;1,3;1,4;2,3;2,4;3,4];
len_multiComp = size(comp_idx,1);
locResp_abs_spatialD = cell(lenM,len_absD); %modality x abs spatialD
[meanLocResp_abs_spatialD, SEM_locResp_abs_spatialD] = deal(NaN(lenM,len_absD));
for i = 1:len_absD %4 abs spatial discrepancy
    for j = 1:lenM
        locResp_abs_spatialD_subj = eval(['locResp',modality{j},...
            '_abs_AVd_subj(i:len_absD:end,:);']);
        locResp_abs_spatialD{j,i} = mean(locResp_abs_spatialD_subj);
        meanLocResp_abs_spatialD(j,i) = mean(locResp_abs_spatialD{j,i});
        SEM_locResp_abs_spatialD(j,i) = std(locResp_abs_spatialD{j,i})/sqrt(length(subjNs));
    end
end

% save it as excel file so we can also do pairwise t-tests in R
locResp_abs_spatialD_temp   = cell2mat(locResp_abs_spatialD);
locResp_abs_spatialD_temp_A = locResp_abs_spatialD_temp(1,:);
locResp_abs_spatialD_temp_V = locResp_abs_spatialD_temp(2,:);
locResp_abs_spatialD_A_reshape = reshape(locResp_abs_spatialD_temp_A,...
    [length(subjNs), len_absD])';
locResp_abs_spatialD_V_reshape = reshape(locResp_abs_spatialD_temp_V,...
    [length(subjNs), len_absD])';
locResp_abs_spatialD_A_mat = locResp_abs_spatialD_A_reshape(:);
locResp_abs_spatialD_V_mat = locResp_abs_spatialD_V_reshape(:);

abs_spatialD_mat = repmat(unique_absD,[length(subjNs), 1]);
subjI_mat        = cell(len_absD, length(subjNs));
for i = 1:length(subjNs)
    for j = 1:len_absD
        subjI_mat{j,i} = subjIs{i};
    end
end
subjI_mat = subjI_mat(:);

Table_VE_absSpatialD = table(locResp_abs_spatialD_A_mat,...
    locResp_abs_spatialD_V_mat, abs_spatialD_mat, subjI_mat,...
    'VariableNames', {'VE_Ashift_absD','VE_Vshift_absD','SpatialD_abs','SubjI'});
writetable(Table_VE_absSpatialD,'Table_VE_absSpatialD.xlsx','Sheet','MyNewSheet',...
    'WriteVariableNames',true);
writetable(Table_VE_absSpatialD,'Table_VE_absSpatialD.txt');

