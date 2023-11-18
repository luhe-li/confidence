%% load data from unimodal localization task
clear all; close all; clc; %parpool(3); %3 cores/12 hours
subjN = 3;
subjI = 'PW';
% addpath(genpath('/home/fh862/bads-master')); %HPC
% addpath(genpath('/home/fh862/Proj_updatePrior_CausalInference')); %HPC
% addpath(genpath(['/home/fh862/Proj_updatePrior_',subjI,'/data_',subjI])); %HPC
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                        'Unimodal localization v3/Data/', subjI]));
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
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                        'Pointing task/Data/', subjI]));
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
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'Matching task v3/Data/',subjI]));
E = load(['AV_alignment_sub', num2str(subjN), '_dataSummary.mat'],...
    'AV_alignment_data');
data.matching = E.AV_alignment_data{end}.sorted;

%% load data from bimodal localization task
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'Bimodal localization v2/Data/',subjI]));
F = load(['BimodalLocalization_sub', num2str(subjN), '_dataSummary.mat'],...
    'BimodalLocalization_data');
%size: 2 (adaptation cond) x 2 (phase) x 4 (A loc) x 4 (V loc)
data.bimodal_unity          = F.BimodalLocalization_data{2};
data.bimodal_unity_prob     = F.BimodalLocalization_data{end-1};
data.bimodal_locResp        = F.BimodalLocalization_data{end};
data.numUnityTrialsPerLoc   = 20;

%% Fit the full model jointly to all the data
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'ModelFitting/Fits']));
Mcmp                 = load('Results_modelComparison_overall.mat', 'ModelComparison');
sub_idx              = find(subjN==Mcmp.ModelComparison{1});
bestM                = Mcmp.ModelComparison{end}{sub_idx};
%e.g., 'SW-samePC1cond-strategyMAP_MA_strategyUnity_posteriorC1'
seg_idx              = find(bestM == '_'); 
model.strategy_MAP   = bestM((seg_idx(1)+1):(seg_idx(2)-1)); %'MS' or 'MA'
model.strategy_unity = bestM((seg_idx(3)+1):end); %measurements, MAP, posteriorC1

%% Fit the full model jointly to all the data
addpath(genpath(['/Users/hff/Desktop/NYU/Project1/Experiment ',...
                        'code_project1/','ModelFitting/bads-master']));
model.numBins_A      = 100;
model.numBins_V      = 100;
model.numSD          = 5;
model.cond           = {'cong','incong'};
model.phase          = {'pre', 'post'};
model.modality       = ['A','V'];

%free parameters for the full model
%--------------------------------------------------------------------------
%p(1)-p(2)  : a_A, b_A, terms that define the relative auditory biases
%p(3)-p(4)  : sigma_A, sigma_V for unimodal trials
%p(5)-p(7)  : sigma_AV_A_pre, sigma_AV_A_post_cong, sigma_AV_A_post_incong,
%p(8)-p(10) : sigma_AV_V_pre, sigma_AV_V_post_cong, sigma_AV_V_post_incong,
%p(11)      : episolon (only applicable for some strategy)
%p(12)-p(13): pC1_cong, pC1_incong, pCommon for the pre- and the
%               post-association phase
%p(14)-p(16): lapse rate for the matching task and the unity judgments for
%               the congruent and incongruent condition
%--------------------------------------------------------------------------
nLogL     = @(p) nLL_overall_samePC1cond_diffSigma(p(1),p(2),p(3),p(4),p(5),p(6),p(7),...
                p(8),p(9),p(10),p(11),p(12),p(13),p(14),p(15),p(16),data,model);
model.lb  = [0.5,-20, 1e-1, 1e-2,1e-1,1e-1,1e-1, 1e-2, 1e-2, 1e-2, 0.1, 1e-2, 1e-2, 1e-4, 1e-4, 1e-4];
model.ub  = [  3, 10,   15,    5,  15,  15,  15,    5,    5,    5,  20, 0.99, 0.99, 0.06, 0.06, 0.06];
model.plb = [1.5, -9,    1,  0.1,   5,   5,   5,    1,    1,    1,   5,  0.1,  0.1, 1e-2, 1e-2, 1e-2];
model.pub = [  2,  0,    5,    2,  10,  10,  10,    3,    3,    3,  10,  0.5, 0.5, 0.03, 0.03, 0.03];
if strcmp(model.strategy_unity,'posteriorC1')
    model.lb(11) = 0; model.ub(11) = 0; model.plb(11) = 0; model.pub(11) = 0;
end
model.numRuns     = 100;
model.numSect_bds = 4;
model.inits       = getInit(model.plb, model.pub, model.numSect_bds, model.numRuns);
%initialize matrices that store the best-fit estimates and minimum nLL.
P       = NaN(model.numRuns, length(model.lb));
minNLL  = NaN(1, model.numRuns);
for i = 1:model.numRuns
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
    '_strategyUnity_', model.strategy_unity, '_diffSigma_sub', num2str(subjN), '_',...
    datestr(datetime('now')), '.mat'],'ModelFitting');
    

