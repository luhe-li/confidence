%% load data from unimodal localization task
clear all; close all; clc; parpool(2); %2 cores, 2 hours
subjN = 15;
subjI = 'AD';
addpath(genpath('/home/fh862/bads-master')); %HPC
addpath(genpath('/home/fh862/Proj_updatePrior_CausalInference')); %HPC
addpath(genpath(['/home/fh862/Proj_updatePrior_',subjI,'/data_',subjI])); %HPC
% addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
%                         'Unimodal localization v3/Data/', subjI]));
C                    = load(['Unimodal_localization_sub', num2str(subjN), '.mat'],...
                        'Unimodal_localization_data');
uniLoc_A             = C.Unimodal_localization_data{end}.data(1:2,:); 
                        %ignore response time (stored in the 3rd row)
data.s_A             = unique(uniLoc_A(1,:));
uniLoc_V             = C.Unimodal_localization_data{end-1}.data(1:2:end,:);
                        %ignore V locations in cm (stored in the 2nd row)
                        %and response time (stored in the 4th row)
data.s_V             = unique(uniLoc_V(1,:));
data.unimodal(:,:,1) = uniLoc_A; 
data.unimodal(:,:,2) = uniLoc_V;

%% load data from the pointing practice
% addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
%                         'Pointing task/Data/', subjI]));
D                    = load(['PointingTest_sub', num2str(subjN), '.mat'],...
                        'PointingTest_data');
cursorLoc            = D.PointingTest_data{end}(1,:);
cursorLoc_unique     = unique(cursorLoc);
cursor_locResp       = D.PointingTest_data{end}(2,:);
%remove outliers
locError             = cursor_locResp - cursorLoc;
SD                   = std(locError);
bool_nonoutliers     = (locError./SD < 3 & locError./SD > -3);
cursor_locResp_rm    = locError(bool_nonoutliers);
%We use MLE to calculate sigma_r, which is sqrt(SSE/N)
%We can't just use function std.m because the denominator is N-1
data.sigma_r         = sqrt(sum((cursor_locResp_rm - mean(cursor_locResp_rm)).^2)/...
                        length(cursor_locResp_rm));

%% load data from the matching task
% addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
%     'Matching task v3/Data/',subjI]));
E = load(['AV_alignment_sub', num2str(subjN), '_dataSummary.mat'],...
    'AV_alignment_data');
data.matching = E.AV_alignment_data{end}.sorted;

%% load data from bimodal localization task
% addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
%     'Bimodal localization v2/Data/',subjI]));
F = load(['BimodalLocalization_sub', num2str(subjN), '_dataSummary.mat'],...
    'BimodalLocalization_data');
%size: 2 (adaptation cond) x 2 (phase) x 4 (A loc) x 4 (V loc)
data.bimodal_unity          = F.BimodalLocalization_data{2};
data.bimodal_unity_prob     = F.BimodalLocalization_data{end-1};
data.bimodal_locResp        = F.BimodalLocalization_data{end};
data.numUnityTrialsPerLoc   = 20;

addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'Adaptation v2/Data/',subjI]));
try
    G  = load(['Adaptation_congruent_sub',num2str(subjN),'_session2.mat'],...
        'Adaptation_data');
    H = load(['Adaptation_incongruent_sub', num2str(subjN),'_session1.mat'],...
        'Adaptation_incongruent_data');
catch
    G  = load(['Adaptation_congruent_sub',num2str(subjN),'_session1.mat'],...
        'Adaptation_data');
    H = load(['Adaptation_incongruent_sub', num2str(subjN),'_session2.mat'],...
        'Adaptation_incongruent_data');
end
%congruent condition
data.learning_sT      = 110:160; %the last 1/3 of the trials during the learning phase
data.cong_Vlocs       = G.Adaptation_data{end-2}.arrangedLocs_deg(data.learning_sT);
data.cong_Alocs       = G.Adaptation_data{end-1}.arrangedLocs_deg(data.learning_sT);
data.cong_locModality = G.Adaptation_data{end-4}.localize_modality(data.learning_sT);
data.cong_locResp     = G.Adaptation_data{end}.localization(data.learning_sT);
data.cong_unity       = G.Adaptation_data{end}.unity(data.learning_sT);

%incongruent condition
data.incong_Vlocs       = H.Adaptation_incongruent_data{end-2}.arrangedLocs_deg(data.learning_sT);
data.incong_Alocs       = H.Adaptation_incongruent_data{end-1}.arrangedLocs_deg(data.learning_sT);
data.incong_locModality = H.Adaptation_incongruent_data{end-4}.localize_modality(data.learning_sT);
data.incong_locResp     = H.Adaptation_incongruent_data{end}.localization(data.learning_sT);
data.incong_unity       = H.Adaptation_incongruent_data{end}.unity(data.learning_sT);

%% Fit the full model jointly to all the data
% addpath(genpath(['/Users/hff/Desktop/NYU/Project1/Experiment code_project1/',...
%                         'ModelFitting/bads-master']));
model.numBins_A      = 100;
model.numBins_V      = 100;
model.numSD          = 5;
model.cond           = {'cong','incong'};
model.phase          = {'pre', 'post'};
model.modality       = ['A','V'];
model.strategy_MAP   = 'MS'; %'MS' or 'MA'
model.strategy_unity = 'measurements'; %measurements, MAP, posteriorC1

%free parameters for the full model
%--------------------------------------------------------------------------
%p(1)-p(2)  : a_A, b_A, terms that define the relative auditory biases
%p(3)-p(4)  : sigma_A, sigma_V for unimodal trials
%p(5)-p(6)  : sigma_AV_A, sigma_AV_V for bimodal trials
%p(7)       : episolon (only applicable for some strategy)
%p(8)-p(9)  : pC1_pre_cong, pC1_post_cong, pCommon for the pre- and the
%               post-association phase
%p(10)-p(11): pC1_pre_incong, pC1_post_incong, pCommon for the pre- and the
%               post-dissociation phase
%p(12)-p(14): lapse rate for the matching task and the unity judgments for
%               the congruent and incongruent condition
%--------------------------------------------------------------------------
nLogL     = @(p) nLL_overall_wLearningPhase(p(1),p(2),p(3),p(4),p(5),p(6),...
            p(7),p(8),p(9),p(10),p(11),p(12),p(13),p(14),data,model);
model.lb  = [0.5,-20, 1e-1, 1e-2,1e-1, 1e-2, 0.1, 1e-2, 1e-2,1e-2, 1e-2, 1e-4, 1e-4, 1e-4];
model.ub  = [3.5, 10,   15,    5,  20,   10,  30, 0.99, 0.99,0.99, 0.99, 0.06, 0.06, 0.06];
model.plb = [1.5, -9,    1,  0.1,   5,    2,   5,  0.1,  0.1, 0.1,  0.1, 1e-2, 1e-2, 1e-2];
model.pub = [  2,  0,    5,    2,  10,    5,  10,  0.5,  0.5, 0.5,  0.5, 0.03, 0.03, 0.03];
if strcmp(model.strategy_unity,'posteriorC1')
    model.lb(7) = 0; model.ub(7) = 0; model.plb(7) = 0; model.pub(7) = 0;
end
model.numRuns     = 20;
model.numSect_bds = 4;
model.inits       = getInit(model.plb, model.pub, model.numSect_bds, model.numRuns);
%initialize matrices that store the best-fit estimates and minimum nLL.
P       = NaN(model.numRuns, length(model.lb));
minNLL  = NaN(1, model.numRuns);
parfor i = 1:model.numRuns
    disp(i)
    try
        [P(i,:),minNLL(i)] = bads(nLogL, model.inits(i,:), model.lb, model.ub,...
                                model.plb, model.pub);
        disp(P(i,:));
    catch
        disp('undefined at the initial points.');
    end
end
%find the best P with the minimum nLL
fits.P                   = P;
fits.minNLL              = minNLL;
[fits.minNLL_f, min_idx] = min(fits.minNLL); %final
fits.P_f                 = P(min_idx,:);

%% save the data
ModelFitting = {subjN, subjI, data, model, fits};
save(['ModelFitting_updatePrior_fullModel_strategyMAP_', model.strategy_MAP, ...
    '_strategyUnity_', model.strategy_unity, '_wLearningPhase_sub', ...
    num2str(subjN), '_', datestr(datetime('now')), '.mat'],'ModelFitting');
    

