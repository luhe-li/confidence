%% load data from unimodal localization task
clear all; close all; clc; parpool(10);
subjN = 24;    %   3,    6,    9,   11,   12,   19,   20,   22,   24
subjI = 'MD';  %'PW', 'YZ', 'ZZ', 'BB', 'ZY', 'RE', 'MM', 'CS', 'MD'
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
data.pointing        = D.PointingTest_data{end}(1:2,:);
cursorLoc            = data.pointing(1,:);
cursorLoc_unique     = unique(cursorLoc);
cursor_locResp       = data.pointing(2,:);
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
data.bimodal_unity        = F.BimodalLocalization_data{2};
data.bimodal_unity_prob   = F.BimodalLocalization_data{end-1};
data.bimodal_locResp      = F.BimodalLocalization_data{end};
data.numUnityTrialsPerLoc = 20;

%% Fit the full model jointly to all the data
% addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
%     'ModelFitting/Fits']));
Mcmp                   = load('Results_modelComparison_overall.mat', 'ModelComparison');
sub_idx                = find(subjN==Mcmp.ModelComparison{1});
bestM                  = Mcmp.ModelComparison{end}{sub_idx};
%e.g., 'SW-full-strategyLocResp_MA_strategyUnity_posteriorC1'
seg_idx                = find(bestM == '_'); 
model.strategy_locResp = bestM((seg_idx(1)+1):(seg_idx(2)-1)); %'MS' or 'MA'
model.strategy_unity   = bestM((seg_idx(3)+1):end); %measurements, MAP, posteriorC1
if strcmp(model.strategy_locResp,'MS') == 1 && strcmp(model.strategy_unity,'MAP') == 1
    disp('Error! It''s circular!');
end

%% Fit the full model jointly to all the data
% addpath(genpath(['/Users/hff/Desktop/NYU/Project1/Experiment ',...
%                         'code_project1/','ModelFitting/bads-master']));
model.numBins_A      = 100;
model.numBins_V      = 100;
model.numSD          = 5;
model.cond           = {'cong','incong'};
model.phase          = {'pre', 'post'};
model.modality       = ['A','V'];

%free parameters for the full model
%--------------------------------------------------------------------------
%p(1)-p(2)  : a_A, b_A, terms that define the relative auditory biases
%p(3)       : sigma_P, the spread of the central prior
%p(4)-p(5)  : sigma_A, sigma_V for unimodal trials
%p(6)-p(7)  : sigma_AV_A, sigma_AV_V for bimodal trials
%p(8)       : episolon (only applicable for some strategy)
%p(9)-p(10) : pC1_pre_cong, pC1_post_cong, pCommon for the pre- and the
%               post-association phase
%p(11)-p(12): pC1_pre_incong, pC1_post_incong, pCommon for the pre- and the
%               post-dissociation phase
%p(13)-p(15): lapse rate for the matching task and the unity judgments for
%               the congruent and incongruent condition
%--------------------------------------------------------------------------
model.lb  = [0.5,-20,  1, 1e-1, 1e-3,1e-1, 1e-3, 0.1, 1e-2, 1e-2,1e-2, 1e-2, 1e-4, 1e-4, 1e-4];
model.ub  = [3.5,  5, 30,   15,    5,  25,   10,  10, 0.99, 0.99,0.99, 0.99, 0.06, 0.06, 0.06];
model.plb = [  1, -5,  5,    1,  0.1,  10,    2,   1,  0.1,  0.1, 0.1,  0.1, 1e-2, 1e-2, 1e-2];
model.pub = [2.5,  0, 20,    5,    2,  25,    5,   5,  0.5,  0.5, 0.5,  0.5, 0.03, 0.03, 0.03];
if strcmp(model.strategy_unity,'posteriorC1')
    model.lb(8) = 0; model.ub(8) = 0; model.plb(8) = 0; model.pub(8) = 0;
end
model.numRuns     = 20;
model.numSect_bds = 4;
model.inits       = getInit(model.plb, model.pub, model.numSect_bds, model.numRuns);%bootstrap data

%% bootstrap data
model.numBtst  = 120;
D_btst         = btst_allData(data.matching, data.pointing, data.unimodal,...
                 data.bimodal_unity, data.bimodal_locResp, model.numBtst);
fits           = cell(1,model.numBtst);

parfor n = 1:model.numBtst
    disp(['Btst = ',num2str(n)]);
    %define negative log likelihood function (D is changing for each n)
    nLogL   = @(p) nLL_HighPlasticityLongLasting(p(1),p(2),p(3),p(4),p(5),...
                    p(6),p(7),p(8),p(9),p(10),p(11),p(12),p(13),p(14),...
                    p(15),D_btst{n},model);
    %initialize matrices that store the best-fit estimates and minimum nLL.
    P       = NaN(model.numRuns, length(model.lb));
    minNLL  = NaN(1,model.numRuns);
    for i = 1:model.numRuns
        disp(i)
        try
            [P(i,:),minNLL(i)] = bads(nLogL, model.inits(i,:), model.lb,...
                model.ub, model.plb, model.pub);
            disp(squeeze(P(i,:)));
        catch
            disp('undefined at the initial points.');
        end
    end
    
    %find the best P with the minimum nLL
    [minNLL_f, min_idx] = min(minNLL); %final
    P_f                 = P(min_idx,:);
    fits{n} = {P, minNLL, minNLL_f, P_f};
end

%% save the data
ModelFitting = {subjN, subjI, data, D_btst, model, fits};
save(['ModelFitting_updatePrior_HighPlasticityLongLasting_strategyLocResp_',...
    model.strategy_locResp, '_strategyUnity_', model.strategy_unity, ...
    '_btst_sub', num2str(subjN), '_', datestr(datetime('now')), '.mat'],...
    'ModelFitting');
    

