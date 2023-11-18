function [p, nLL] = plotModelFits(modelTypes, stratG, bool_plot, bool_save)
%---------------------------------------------------------------------------
%subjN, subjI: identity of 1 subject; we can't pass multiple subjects once
%modelTypes  : 'fullModel', 'samePC1_pre', 'samePC1_cond','nullModel'
%bool_data   : [1,0] unity judgments were included in the model fitting
%              [0,1] localization responses were included in the model fitting
%              [1,1] both datasets were included
%bool_plot   : we can plot (1) the matching task
%                          (2) the proportion of reporting C=1 with predictions
%                          (3) the localization responses with predictions
%bool_plot and bool_save are a vector of booleans.
%---------------------------------------------------------------------------
if nargin < 3; bool_save = [0,0,0]; 
elseif nargin < 2; bool_plot = [0,1,0]; bool_save = [0,0,0]; end

%MT: models of common-cause prior; DS: models of decision strategies
global lenMT lenDS 

%first load data, model info, best-fit p, minimum LL, sigma_r, s_A and s_V
[data, modelInfo_cell, p, nLL] = loadModelFits(modelTypes, stratG);

for i = 1:lenMT
    for j = 1:lenDS
        modelInfo = modelInfo_cell{i,j};      
        %calculate model predictions for each model
        %R1: mu_shat_uni, sigma_shat_wN
        %R2: pC1_given_sAsV, p_MAP, p_MAP_edges, counts
        [R1, R2, pC1] = getModelPredictions(p{i,j}, modelTypes{i}, data, modelInfo);
        
        %plot unity judgment
        if bool_plot(1) == 1 
            %compute the proportion of reporting C=1 predicted by the model
            %and the proportion of reporting C=1 based on parametric
            %bootstrap
            [propC1_modelPred, propC1_param_btst_68CI, propC1_abs_modelPred,...
                propC1_abs_param_btst_68CI] = predictPropC1(R2);
            plotUnity(data.s_V, modelInfo, propC1_modelPred, propC1_param_btst_68CI,...
                propC1_abs_modelPred, propC1_abs_param_btst_68CI, pC1,...
                modelTypes{i}, stratG{j}, bool_save(1));
        end
        
        %plot indivisual localization responses
        if bool_plot(3) == 1 || bool_plot(2) == 1
            [s,p_MAP] = plotlocResp(data.s_V, R2, data, pC1, modelInfo,...
                modelTypes{i},stratG{j}, bool_plot(2), bool_save(2)); 
        end
        
        %plot the ventriloquism effects
        if bool_plot(3) == 1           
            %compute the ventriloquism effects predicted by the model
            %those based on parametric bootstrap
            [VE_modelPred, VE_param_btst_95CI, VE_abs_modelPred,...
                VE_abs_param_btst_95CI] = predictVE(s, p_MAP, R1);
            plotVE(data.s_V, VE_modelPred, VE_param_btst_95CI, VE_abs_modelPred, ...
                VE_abs_param_btst_95CI, modelTypes{i}, bool_save(3)); 
        end
    end
end

%==========================================================================
%                        MODEL PREDICTIONS
%==========================================================================
function [data, MInfo, bestP, minNLL] = loadModelFits(M, stratG)
%output--------------------------------------------------------------------
%data: data from all the three preparatory and the main experiment
%MInfo: details for how each model was fitted to the data (e.g.,
%       initializations, lb, ub and etc)
%bestP: estimated parameter given the best-fitting model
%minNLL: the minimum negative log likelihood

%input---------------------------------------------------------------------
%M: best-fitting model of 
%stratG: best-fitting strategy

global sN sI lenMT lenDS
%add path and initialize a cell and a matrix
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'ModelFitting/Fits/', sI]))

[bestP, MInfo] = deal(cell(lenMT,lenDS)); 
minNLL         = NaN(lenMT,lenDS);
for i = 1:lenMT
    for j = 1:lenDS
        %first load the file
        C = load(['ModelFitting_updatePrior_',M{i},'_',stratG{j},'_sub',...
            num2str(sN),'.mat'], 'ModelFitting');
        %load the data
        if i == 1; data = C.ModelFitting{3}; end
        %load the best-fit parameters and minNLL for each model
        bestP{i,j}  = C.ModelFitting{end}.P_f;
        minNLL(i,j) = C.ModelFitting{end}.minNLL_f;
        MInfo{i,j}  = C.ModelFitting{4};
    end
end

function [R1, R2, pC1] = getModelPredictions(p, M, data, modelInfo)
%output--------------------------------------------------------------------
%R1: mu_shat_uni, sigma_shat_wN
%R2: model predicted distributions for loc responses in bimodal
%pC1: the common-cause prior

%input---------------------------------------------------------------------
%p: estimated parameters corresponding to minimum nLL
%M: best-fitting model
%data: behavioral data
%ModelInfo: how the model is fitted, e.g., number of fits, initializations

%calculate model prediction by passing the best-fit parameters to nLL funciton
switch M
    case 'fullModel' 
        [~, R1, R2] = nLL_overall(p(1),p(2),p(3),p(4),p(5),p(6),p(7),...
            p(8),p(9),p(10),p(11),p(12),p(13),p(14),data,modelInfo);
        pC1 = round([p(8), p(9); p(10), p(11)],2);
    case 'samePC1pre'
        [~, R1, R2] = nLL_overall_samePC1pre(p(1),p(2),p(3),p(4),p(5),...
            p(6),p(7),p(8),p(9),p(10),p(11),p(12),p(13),data,modelInfo);
        pC1 = round([p(8), p(9); p(8), p(10)],2);
    case 'samePC1cond'
        [~, R1, R2] = nLL_overall_samePC1cond(p(1),p(2),p(3),p(4),p(5),...
            p(6),p(7),p(8),p(9),p(10),p(11),p(12),data,modelInfo);
        pC1 = round([p(8), p(8); p(9), p(9)],2);      
    case 'null'
        [~, R1, R2] = nLL_null(p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),...
            p(9),p(10),p(11),data,modelInfo);
        pC1 = round([p(8), p(8); p(8), p(8)],2);      
end

function [propC1_pred, propC1_param_btst_68CI, propC1_abs_pred, ...
    propC1_abs_param_btst_68CI] = predictPropC1(R2)
%output--------------------------------------------------------------------
%propC1_pred: the proportion of reporting C=1 as a function of spatial
%           discrepancy predicted by a given model
%propC1_param_btst_68CI: the 68% confidence interval (1 SD) on the 
%           proportion of reporting C=1 as a function of spatial 
%           discrepancy by parametric bootstrap
%propC1_abs_pred: the proportion of reporting C=1 as a function of absolute 
%           spatial discrepancy predicted by a given model
%propC1_abs_param_btst_68CI: the 68% confidence interval (1 SD) on the 
%           proportion of reporting C=1 as a function of absolute spatial 
%           discrepancy by parametric bootstrap

%input---------------------------------------------------------------------
%R2: model predicted distributions for loc responses in bimodal

global lenC lenP lenD lenS
%number of trials for each spatial discrepancy (V-A=-24,-16,-8,8,16,24)
numT_AV                    = [20,40,60,80,60,40,20];
%number of parametric simulations
n_param_btst               = 1e3;
%the indices for 16% and 84% percentiles
idx_CI_lb                  = floor(n_param_btst*0.16);
idx_CI_ub                  = ceil(n_param_btst*0.84);
%initialize the proportion of reporting C=1
propC1_pred                = NaN(lenC, lenP, lenD);
propC1_abs_pred            = NaN(lenC, lenP, ceil(lenD/2));
propC1_param_btst_68CI     = NaN(lenC, lenP, lenD, 2); %last dimension:lb, ub
propC1_abs_param_btst_68CI = NaN(lenC, lenP, ceil(lenD/2), 2);

for i = 1:lenC %for each condition
    %model prediction
    [overall_pC1_pre_model, overall_pC1_post_model] = deal(NaN(lenS,lenD));
    for k = 1:lenS 
        %get the indices for where to store the data
        idx_lb = ceil(lenD/2)+1-k;
        idx_ub = lenD + 1 -k;
        %by model predictions
        overall_pC1_pre_model(k,idx_lb:idx_ub)  = squeeze(R2{i,1}.pC1_given_sAsV(k,:))';
        overall_pC1_post_model(k,idx_lb:idx_ub) = squeeze(R2{i,2}.pC1_given_sAsV(k,:))';
    end
    %average across trials with the same spatial discrepancy
    propC1_pred(i,1,:) = nanmean(overall_pC1_pre_model);
    propC1_pred(i,2,:) = nanmean(overall_pC1_post_model);
    
    %parametric bootstrap
    for j = 1:lenP    
        %selecting spatial discrepancy = 0, 8, 16, 24
        propC1_abs_pred(i,j,:)     = propC1_pred(i,j,ceil(lenD/2):end);
        %taking an average of spatial discrepancy of +/-8, +/-16, +/-24
        propC1_abs_pred(i,j,2:end) = propC1_abs_pred(i,j,2:end)./2 + ...
                                     propC1_pred(i,j,floor(lenD/2):-1:1)./2;
        %initialize propC1_param_btst for condition i and phase j
        propC1_param_btst_ij = NaN(lenD, n_param_btst);
        for l = 1:lenD
            %given the model-predicted distribution, draw # samples 
            %(binary: 0 or 1)
            rBinary_param_btst = parametric_simulation([0,1], ...
                [1-propC1_pred(i,j,l),propC1_pred(i,j,l)],  numT_AV(l),...
                n_param_btst);
            %compute the mean across trials with the same absolute spatial
            %discrepancy
            propC1_param_btst_ij(l,:)       = mean(rBinary_param_btst, 2);
            %sort all 1000 mean  in an ascending order
            propC1_param_btst_sort          = sort(propC1_param_btst_ij(l,:));
            %get the confidence interval
            propC1_param_btst_68CI(i,j,l,1) = propC1_param_btst_sort(idx_CI_lb);
            propC1_param_btst_68CI(i,j,l,2) = propC1_param_btst_sort(idx_CI_ub);
        end
        %for computing propC1 as a function of absolute spatial discrepancy
        propC1_abs_param_btst = propC1_param_btst_ij(ceil(lenD/2):end,:);
        propC1_abs_param_btst(2:end,:) = propC1_abs_param_btst(2:end,:)./2 + ...
                                         propC1_param_btst_ij(floor(lenD/2):-1:1,:)./2;
        propC1_abs_param_btst_sort          = sort(propC1_abs_param_btst,2);
        propC1_abs_param_btst_68CI(i,j,:,1) = propC1_abs_param_btst_sort(:,idx_CI_lb);
        propC1_abs_param_btst_68CI(i,j,:,2) = propC1_abs_param_btst_sort(:,idx_CI_ub);
    end
end

function [VE_modelPred, VE_param_btst_95CI, VE_abs_modelPred, ...
    VE_abs_param_btst_95CI] = predictVE(s, p_MAP, R1)
%output--------------------------------------------------------------------
%VE_modelPred: ventriloquism effects as a function of spatial discrepancy
%               predicted by the best-fitting model
%VE_param_btst_95CI: 95% confidence interval on ventriloquism effects as a 
%               function of spatial discrepancy using parametric bootstrap
%VE_abs_modelPred: ventriloquism effects as a function of absolute spatial 
%               discrepancy predicted by the best-fitting model
%VE_abs_param_btst_95CI: 95% confidence interval on ventriloquism effects 
%               as a function of absolute spatial discrepancy using 
%               parametric bootstrap

%input---------------------------------------------------------------------
%s: stimulus location (a wide range with a fine bin size)
%p_MAP: the model-predicted probability 
%       size(p_MAP): lenC x lenP (length(s_A) x length(s_V) x lenM x length(s))
%R1.mu_shat_uni: model-predicted perceived stimulus location when presented
%       unimodally
%R1.sigma_shat_wN: measurement noise of an auditory stimulus when presented alone

global lenD lenS lenC lenP lenM
%number of localization trials for each audiovisual pair
nT_bimodal_perM        = 10;
%number of localization trials for each unimodal stimulus at each location
nT_unimodal_perM       = 30;
%number of parametric simulations
n_param_btst           = 1e3;
%the indices for 2.5% and 97.5% percentiles
idx_lb                 = floor(n_param_btst*0.025);
idx_ub                 = ceil(n_param_btst*0.975);
%initialization
VE_modelPred           = NaN(lenC, lenP, lenM, lenD);
VE_param_btst_95CI     = NaN(lenC, lenP, lenM, lenD, 2); %last dimension: lb, ub
VE_abs_modelPred       = NaN(lenC, lenP, lenM, ceil(lenD/2));
VE_abs_param_btst_95CI = NaN(lenC, lenP, lenM, ceil(lenD/2), 2);
r_uni_param_btst       = NaN(lenS, lenM, n_param_btst);

%parametric bootstrapping unimodal trials
for k = 1:lenS
    for l = 1:lenM
        if l == 1; mu_uni = R1.mu_shat_A_uni(k); sig_uni = R1.sigma_shat_A_wN;
        else; mu_uni = R1.mu_shat_V_uni(k); sig_uni = R1.sigma_shat_V_wN;end
        %analytical solution
        uni_dist  = norm_dst(s, mu_uni, sig_uni, 0);
        r_uni_param_btst_kj = parametric_simulation(s,...
            uni_dist, nT_unimodal_perM, n_param_btst);
        %drawing samples using parametric bootstrapping
        r_uni_param_btst(k,l,:) = mean(r_uni_param_btst_kj,2);
    end
end

%parametric bootstrapping bimodal trials
for i = 1:lenC
    for j = 1:lenP
        %select p_MAP for condition i and phase j
        p_MAP_mat = p_MAP{i,j};
        %initialization
        VE_ij            = NaN(lenS,lenM,lenD);
        VE_param_btst_ij = NaN(lenS, lenM, lenD, n_param_btst);
        
        for k = 1:lenS %A loc 
            for l = 1:lenS %V loc
                idx = ceil(lenD/2)-k+l;
                for m = 1:lenM
                    if m == 1; mu_uni = R1.mu_shat_A_uni(k); 
                    else; mu_uni = R1.mu_shat_V_uni(l); end
                    
                    %distribution for localization responses (bimodal)
                    p_MAP_norm           = squeeze(p_MAP_mat(k,l,m,:))';
                    p_MAP_norm           = p_MAP_norm./sum(p_MAP_norm);
                    
                    %VE model prediction (size: 4 locs x 7 spatial discrepancies)
                    % = E[responded stimulus location] in bimodal trials
                    %   - the perceived location in unimodal trials
                    VE_ij(k,m,idx)       = p_MAP_norm*s' - mu_uni;
                    
                    %simulate bimodal auditory localization responses
                    r_bimodal_param_btst = parametric_simulation(s, p_MAP_norm,...
                                           nT_bimodal_perM, n_param_btst);
                    %VE by parametric bootstrapping 
                    if m == 1
                        VE_param_btst_ij(k,m,idx,:) = mean(r_bimodal_param_btst - ...
                            squeeze(r_uni_param_btst(k,m,:)),2);
                    else
                        VE_param_btst_ij(k,m,idx,:) = mean(r_bimodal_param_btst - ...
                            squeeze(r_uni_param_btst(l,m,:)),2);
                    end
                end 
            end
        end
        
        %model prediction & parametric bootstrapping
        for m = 1:lenM
            %here we need to correct the sign:
            %A: positive for spatialD (V-A) = 0,8,16,24; 
            %   negative for spatialD (V-A) = -8,-16,24;
            %V: negative for spatialD (V-A) = 0,8,16,24; 
            %   positive for spatialD (V-A) = -8,-16,24;
            if m == 1; corr_sign = [1,-1]; else; corr_sign = [-1,1]; end
            
            %compute ventriloquism effects based on the model prediction
            %(deterministic)
            VE_modelPred(i,j,m,:)         = nanmean(VE_ij(:,m,:));
            
            %compute ventriloquism effects based as a function of absolute
            %spatial discrepancy based on the model prediction
            VE_abs_modelPred(i,j,m,:)     = VE_modelPred(i,j,m,ceil(lenD/2):end); 
            VE_abs_modelPred(i,j,m,2:end) = corr_sign(1).*VE_abs_modelPred(i,j,m,2:end)./2 +...
                corr_sign(2).*VE_modelPred(i,j,m,floor(lenD/2):-1:1)./2;
                                        
            %compute ventriloquism effects using parametric simulation
            VE_param_btst_nanmean     = squeeze(nanmean(VE_param_btst_ij(:,m,:,:)));
            VE_param_btst_sort        = sort(VE_param_btst_nanmean,2);
            
            %compute ventriloquism effects as a function of absolute spatial discrepancy
            VE_abs_param_btst_nanmean = VE_param_btst_nanmean(ceil(lenD/2):end, :);
            VE_abs_param_btst_nanmean(2:end,:) = corr_sign(1).*...
                VE_abs_param_btst_nanmean(2:end,:)./2 +corr_sign(2).*...
                VE_param_btst_nanmean(floor(lenD/2):-1:1, :)./2;
            VE_abs_param_btst_sort    = sort(VE_abs_param_btst_nanmean,2);
            
            %lower bound of the 95% of the parametric simulation
            VE_param_btst_95CI(i,j,m,:,1)     = VE_param_btst_sort(:, idx_lb);
            VE_abs_param_btst_95CI(i,j,m,:,1) = VE_abs_param_btst_sort(:, idx_lb);
            %upper bound of the 95% of the parametric simulation
            VE_param_btst_95CI(i,j,m,:,2)     = VE_param_btst_sort(:, idx_ub); 
            VE_abs_param_btst_95CI(i,j,m,:,2) = VE_abs_param_btst_sort(:, idx_ub);
        end
    end
end


%==========================================================================
%                         PLOTTING FUNCTIONS
%==========================================================================
function plotVE(s_V, VE_modelPred, VE_95CI_paramSim,VE_abs_modelPred,...
    VE_abs_param_btst_95CI, M, bool_save, bool_plt_95CI_paramSim)
%if not specified, we will plot error bars on the model predictions
if nargin < 8; bool_plt_95CI_paramSim = 1;end
%global variables we need to use
global sI sN lenM lenC lenP lenD; 
%load data from the pre- and the post-learning phase
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                'Bimodal localization v2/Data/',sI]));
C              = load(['Summary_VE_sub', num2str(sN),'.mat'], 'Summary_VE');
%mean ventriloquism effect as a function of spatial discrepancy
meanVE         = C.Summary_VE{2};
%mean ventriloquism effect as a function of absolute spatial discrepancy
meanVE_abs     = C.Summary_VE{5};
%standard error of the mean
seLocR         = C.Summary_VE{4};
seLocR_abs     = C.Summary_VE{7};
%number of trials for each AV discrepancy
numT_AV        = [20,40,60,80,60,40,20];
%number of trials for each absolute spatial discrepancy [80, 120, 80, 40]
numT_abs_AV    = numT_AV(ceil(lenD/2):end) + [0,fliplr(numT_AV(1:floor(lenD/2)))]; 
%compute the level of absolute spatial discrepancy
AV_combs       = combvec(s_V, s_V);
AV_offset      = sort(unique(AV_combs(2,:) - AV_combs(1,:)));
AV_abs_offset  = unique(abs(AV_offset));
%4 levels in total
lenD_abs       = length(AV_abs_offset);

%define variables for the plots
modality       = {'A','V'};
cond           = {'Congruent', 'Incongruent'};
y_lbl          = {'Auditory localization shifts (deg)', ...
                  'Visual localization shifts (deg)'};
lgd_cond_phase = {'Pre-learning', 'Post-learning'};
x_bds_VE       = [AV_offset(1)-3, AV_offset(end)+3];
x_ticks_VE     = [x_bds_VE(1), AV_offset, x_bds_VE(end)];
y_bds_VE       = [AV_offset(1)-2, AV_offset(end)+2];
y_ticks_VE     = AV_offset(find(y_bds_VE(1) < AV_offset,1):...
                           find(y_bds_VE(end) > AV_offset,1,'last'));
cb             = [65,105,225]./255; cr = [0.85, 0.23, 0.28];
cMap_VE        = {min(cb.*2,1),cb; [247,191,190]./255, cr};
jitter         = [0.8,-0.8]; %or [-0.8,0.8];
lw             = 3; %lineWidth
fs_lbls        = 25;
fs_lgds        = 20;

%plot VE as a function of spatial discrepancy
for i = 1:lenM
    figure
    for j = 1:lenC
        subplot(1,lenC,j)
        addBackground(x_bds_VE, y_bds_VE, x_ticks_VE, [y_bds_VE(1),y_ticks_VE,y_bds_VE(end)])
        %plot a diagonal line indicating complete capture by the other modality
        if i == 1
            plot([AV_offset(1),AV_offset(end)],[AV_offset(1),AV_offset(end)],...
                'k--','Color',ones(1,3).*0.5,'lineWidth',lw); hold on;
        else
            plot([AV_offset(1),AV_offset(end)],[AV_offset(end),AV_offset(1)],...
                'k--','Color',ones(1,3).*0.5,'lineWidth',lw); hold on;
        end
        %plot a horizontal line indicating no ventriloquism effect
        plot(AV_offset, zeros(1,lenD), 'k:','Color',ones(1,3).*0.5,'lineWidth',lw); 
        for k = 1:lenP
            %plot model prediction
            plot(AV_offset+jitter(k), squeeze(VE_modelPred(j,k,i,:)),'-',...
                'lineWidth',lw,'Color',cMap_VE{i,k});
            %plot confidence interval computed by parametric simulation
            if bool_plt_95CI_paramSim
                if k == 1; alpha_vec = [0.3, 0.6]; else; alpha_vec = [0.15, 0.3]; end
                VE_95CI_paramSim_lb = squeeze(VE_95CI_paramSim(j,k,i,:,1));
                VE_95CI_paramSim_ub = squeeze(VE_95CI_paramSim(j,k,i,:,2));
                patch([AV_offset+jitter(k), fliplr(AV_offset+jitter(k))], ...
                    [VE_95CI_paramSim_lb',flipud(VE_95CI_paramSim_ub)'],...
                    cMap_VE{i,k},'FaceAlpha', alpha_vec(1), 'EdgeColor',...
                    cMap_VE{i,k},'EdgeAlpha', alpha_vec(2),'lineWidth',2); 
            end
        end
        
        for k = 1:lenP    
            %plot behavioral data
            for n = 1:lenD
                plt(k) = errorbar(AV_offset(n)+jitter(k), meanVE(j,k,i,n),...
                    2.*seLocR(j,k,i,n), '-o','MarkerSize',sqrt(numT_AV(n)).*2,...
                    'Color',cMap_VE{i,k},'MarkerFaceColor',cMap_VE{i,k},...
                    'MarkerEdgeColor',cMap_VE{i,k},'lineWidth',lw); 
            end
        end
        text(x_bds_VE(1) + 7, y_bds_VE(1)+3, ['Condition: ', cond{j}],'FontSize',fs_lgds); 
        hold off;box off;
        %add legends
        if j == 1; text(x_bds_VE(1) + 0.2, y_bds_VE(end)-1, sI,'FontSize',fs_lgds); end
        legend([plt(1) plt(2)], lgd_cond_phase,'Location','southeast','FontSize',fs_lgds); 
        legend boxoff; xticks(AV_offset); xlim(x_bds_VE);
        xlabel('Spatial discrepancy (V-A, deg)');
        yticks(y_ticks_VE); ylim(y_bds_VE);
        ylabel(sprintf(y_lbl{i})); set(gca,'FontSize',fs_lbls);
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.75, 0.60]);
    set(gcf,'PaperUnits','centimeters','PaperSize',[50 25]);
    if bool_save == 1; saveas(gcf, ['VE_fits_', M,'_cond_',modality{i}, '_', sI], 'pdf'); end
end

%define variables for the new plots
y_bds_VE_abs   = [-8, AV_abs_offset(end)+2];
y_ticks_VE_abs = AV_abs_offset;
x_bds_VE_abs   = [-3, AV_abs_offset(end)+3];
x_ticks_VE_abs = [-3, AV_abs_offset, AV_abs_offset(end)+3];

%plot VE as a function of absolute spatial discrepancy
for i = 1:lenM
    figure
    for j = 1:lenC
        subplot(1,lenM,j)
        addBackground(x_bds_VE_abs, y_bds_VE_abs, x_ticks_VE_abs, [y_bds_VE_abs(1),...
            y_ticks_VE_abs,y_bds_VE_abs(end)])
        %plot a diagonal line indicating complete capture by the other modality
        if i == 1
            plt3 = plot(AV_abs_offset,AV_abs_offset,'k--',...
                'Color',ones(1,3).*0.5,'lineWidth',lw); hold on;
        else
            plt3 = plot(AV_abs_offset,AV_abs_offset,'k--',...
                'Color',ones(1,3).*0.5,'lineWidth',lw); hold on;
        end
        %plot a horizontal line indicating no ventriloquism effect
        plt4 = plot(AV_abs_offset, zeros(1,lenD_abs), 'k:','Color',...
            ones(1,3).*0.5,'lineWidth',lw); hold on;
        for k = 1:lenP
            %plot model prediction
            plot(AV_abs_offset+jitter(k), squeeze(VE_abs_modelPred(j,k,i,:)),...
                '-', 'lineWidth',lw,'Color',cMap_VE{i,k}); hold on
            %plot confidence interval computed by parametric simulation
            if bool_plt_95CI_paramSim
                if k == 1; alpha_vec = [0.3, 0.3]; else; alpha_vec = [0.15, 0.15]; end
                VE_abs_95CI_paramSim_lb = squeeze(VE_abs_param_btst_95CI(j,k,i,:,1));
                VE_abs_95CI_paramSim_ub = squeeze(VE_abs_param_btst_95CI(j,k,i,:,2));
                patch([AV_abs_offset+jitter(k), fliplr(AV_abs_offset+jitter(k))], ...
                    [VE_abs_95CI_paramSim_lb',flipud(VE_abs_95CI_paramSim_ub)'],...
                    cMap_VE{i,k},'FaceAlpha', alpha_vec(1), 'EdgeColor',...
                    cMap_VE{i,k},'EdgeAlpha', alpha_vec(2),'lineWidth',2); hold on;
            end
        end
        
        for k = 1:lenP
            %plot behavioral data
            for n = 1:lenD_abs
                plt(k) = errorbar(AV_abs_offset(n)+jitter(k),...
                    meanVE_abs(j,k,i,n), 2.*seLocR_abs(j,k,i,n),...
                    '-o','MarkerSize',sqrt(numT_abs_AV(n)).*2, 'Color',cMap_VE{i,k},...
                    'MarkerFaceColor',cMap_VE{i,k}, 'MarkerEdgeColor',...
                    cMap_VE{i,k},'lineWidth',lw,'CapSize',15); hold on;
            end
        end
        text(x_bds_VE_abs(1) + 7, y_bds_VE_abs(1)+3, ['Condition: ', cond{j}],...
            'FontSize',fs_lgds); hold off;box off;
        %add legends
        if j == 1; text(x_bds_VE_abs(1) + 0.2, y_bds_VE_abs(end)-1, sI,'FontSize',fs_lgds); end
        legend([plt(1) plt(2)], lgd_cond_phase,'Location','northeast',... %'Position',lgd_pos(j,:),...
            'FontSize',fs_lgds); legend boxoff;
        xticks(AV_abs_offset); xlim(x_bds_VE_abs);
        xlabel('Absolute spatial discrepancy (deg)');
        yticks(y_ticks_VE_abs); ylim(y_bds_VE_abs);
        ylabel(sprintf(y_lbl{i})); set(gca,'FontSize',fs_lbls);
    end
    %best size
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.65, 0.60]);
    set(gcf,'PaperUnits','centimeters','PaperSize',[50 25]); 
    if bool_save == 1; saveas(gcf, ['VE_abs_spatialD_fits_', M,'_cond_',modality{i}, '_', sI], 'pdf'); end
end

function plotUnity(s_V, model, propC1_modelPred, propC1_param_btst_68CI,...
    propC1_abs_modelPred, propC1_abs_param_btst_68CI, pC1, M, G, bool_save,...
    bool_plt_68CI_paramSim)
%if not specified, we will plot error bars on the model predictions
if nargin < 11; bool_plt_68CI_paramSim = 1;end
%global variables we need to use
global sI sN lenC lenP 
%load data from the pre- and the post-learning phase
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                'Bimodal localization v2/Data/',sI]));
C              = load(['Summary_VE_sub', num2str(sN),'.mat'], 'Summary_VE');
%the proportion of reporting 'common-cause' as a function of spatial D
mean_pC1_lenD        = C.Summary_VE{end-3};
%the proportion of reporting 'common-cause' as a function of absolute spatial D
meanPC1_abs_spatialD = C.Summary_VE{end-1};
%SD = sqrt(p*(1-p)/n)
SD_pC1_lenD          = C.Summary_VE{end-2};
SD_pC1_abs_spatialD  = C.Summary_VE{end};

%number of trials for each AV discrepancy
numT_AV        = [20,40,60,80,60,40,20];
%number of trials for each absolute spatial discrepancy [80, 120, 80, 40]
numT_abs_AV    = numT_AV(4:end) + [0,numT_AV(3:-1:1)];
%compute the level of absolute spatial discrepancy
AV_combs       = combvec(s_V, s_V);
AV_offset      = sort(unique(AV_combs(2,:) - AV_combs(1,:)));
AV_abs_offset  = unique(abs(AV_offset));
%4 levels in total
lenD_abs       = length(AV_abs_offset);

%define variables for the plots
x_bds          = [AV_offset(1) - 2, AV_offset(end) + 2]; y_bds = [-0.05,1.05];
x_ticks        = AV_offset; y_ticks = 0:0.25:1;
x_bds_abs      = [-3, AV_abs_offset(end)+3];
x_ticks_abs    = [-3, AV_abs_offset, AV_abs_offset(end)+3];
colorMap       = [0.65, 0.65,0.65;0.1,0.1,0.1];%[135,206,235; 65,105,225]./255;
jitter         = [0.8,-0.8]; %or [-0.8,0.8];
lw             = 3; %lineWidth
fs_lbls        = 25;
fs_lgds        = 20;
lgd_cond_phase = {'Pre-learning', 'Post-learning'};

%plot the unity judgment as a function of spatial discrepancy
figure
for i = 1:lenC
    subplot(1,2,i)
    addBackground(x_bds, y_bds, [x_bds(1),x_ticks, x_bds(end)], [y_bds(1), y_ticks, y_bds(end)])
    %plot the model predictions
    for j = 1:lenP
        plot(AV_offset+jitter(j), squeeze(propC1_modelPred(i,j,:)),'Color',colorMap(j,:),...
            'lineWidth', lw, 'lineStyle','-'); hold on;
        if bool_plt_68CI_paramSim
            propC1_param_btst_68CI_lb = squeeze(propC1_param_btst_68CI(i,j,:,1));
            propC1_param_btst_68CI_ub = squeeze(propC1_param_btst_68CI(i,j,:,2));
            patch([AV_offset+jitter(j), fliplr(AV_offset+jitter(j))], ...
                [propC1_param_btst_68CI_lb',flipud(propC1_param_btst_68CI_ub)'],...
                colorMap(j,:),'FaceAlpha', 0.2, 'EdgeColor',...
                colorMap(j,:),'EdgeAlpha', 0.2,'lineWidth',2); hold on;
        end
    end
    %behavioral data
    for m = 1:length(AV_offset)
        plt_1 = errorbar(AV_offset(m)+jitter(1), mean_pC1_lenD(i,1,m),...
            SD_pC1_lenD(i,1,m),'-o','MarkerSize',sqrt(numT_AV(m)).*2,'Color',...
            colorMap(1,:),'MarkerFaceColor',colorMap(1,:),'MarkerEdgeColor',...
            colorMap(1,:),'lineWidth',lw); hold on
        plt_2 = errorbar(AV_offset(m)+jitter(2),mean_pC1_lenD(i,2,m),...
            SD_pC1_lenD(i,2,m),'-o','MarkerSize',sqrt(numT_AV(m)).*2,'Color',...
            colorMap(2,:),'MarkerFaceColor',colorMap(2,:),'MarkerEdgeColor',...
            colorMap(2,:),'lineWidth',lw); hold on
    end
    
    %add legends
    if i == 1; text(x_bds(1) + 0.2, 1.02, sI,'FontSize',fs_lgds); end
    if i == 1; legend([plt_1 plt_2], lgd_cond_phase,...
            'Position',[0.3 0.15 0.05 0.1]); legend boxoff; end 
    xticks(x_ticks); xlim(x_bds); xlabel('Spatial discrepancy (V - A, deg)'); 
    ylabel(sprintf(['The proportion of reporting a common cause'])); 
    yticks(y_ticks);ylim(y_bds);
    set(gca,'FontSize',fs_lbls);
    %add title
    title([model.cond{i}, 'ruent: p_{C=1,pre} = ', num2str(pC1(i,1)),...
        ', p_{C=1,post} = ', num2str(pC1(i,2))],'FontSize',fs_lgds);
end
M(M=='_')='-'; G(G=='_')='-';
sgtitle([M,' , ',G], 'FontSize',fs_lgds);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.75, 0.60]);
set(gcf,'PaperUnits','centimeters','PaperSize',[50 25]);
if bool_save == 1; saveFig('unityJdg', M, G); end

%plot the unity judgment as a function of absolute spatial discrepancy
figure
for i = 1:lenC
    subplot(1,2,i)
    addBackground(x_bds_abs, y_bds, x_ticks_abs, [y_bds(1), y_ticks, y_bds(end)])
    for j = 1:lenP
        %plot model predictions
        plot(AV_abs_offset+jitter(j), squeeze(propC1_abs_modelPred(i,j,:)),...
            '-','lineWidth',lw,'Color',colorMap(j,:)); hold on;
        %plot parametric bootstrap
        if bool_plt_68CI_paramSim
            propC1_abs_param_btst_68CI_lb = squeeze(propC1_abs_param_btst_68CI(i,j,:,1));
            propC1_abs_param_btst_68CI_ub = squeeze(propC1_abs_param_btst_68CI(i,j,:,2));
            patch([AV_abs_offset+jitter(j), fliplr(AV_abs_offset+jitter(j))], ...
                [propC1_abs_param_btst_68CI_lb',flipud(propC1_abs_param_btst_68CI_ub)'],...
                colorMap(j,:),'FaceAlpha', 0.2, 'EdgeColor',...
                colorMap(j,:),'EdgeAlpha', 0.2,'lineWidth',2); hold on;
        end
    end
    
    for j = 1:lenP
        %behavioral data
        for m = 1:lenD_abs
            plt(j) = errorbar(AV_abs_offset(m)+jitter(j), ...
                meanPC1_abs_spatialD(i,j,m), SD_pC1_abs_spatialD(i,j,m),'-o',...
                'MarkerSize',sqrt(numT_abs_AV(m)).*2,'Color',colorMap(j,:),...
                'MarkerFaceColor',colorMap(j,:),'MarkerEdgeColor',...
                colorMap(j,:),'lineWidth',lw,'CapSize',15); hold on
        end
    end
    %add legends
    if i == 1;text(x_bds_abs(1) + 0.2, 1.02, sI,'FontSize',fs_lgds); hold on;end
    legend([plt(1) plt(2)], lgd_cond_phase,'Location', 'northeast',...
        'FontSize',fs_lgds); legend boxoff;
    xticks(AV_abs_offset); xlim(x_bds_abs); xlabel('Absolute spatial discrepancy (deg)'); 
    ylabel(sprintf(['The proportion of reporting a common cause'])); 
    yticks(y_ticks);ylim(y_bds);
    set(gca,'FontSize',fs_lbls);
end
M(M=='_')='-'; G(G=='_')='-';
sgtitle([M,' , ',G], 'FontSize',fs_lgds);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.65, 0.60]);
set(gcf,'PaperUnits','centimeters','PaperSize',[50 25]);
if bool_save == 1; saveFig('unityJdg_abs_spatialD', M, G); end

function [x,p_MAP] = plotlocResp(s_V, R2, D, pC1, model, M, G, bool_plot, bool_save)
%global variables we need to use
global lenC lenP lenM lenS
%define variables for the plots
cb              = [65,105,225]./255; 
cr              = [0.85, 0.23, 0.28];
cMAP            = {cb, cr};
locResp         = D.bimodal_locResp;
sigma_r         = D.sigma_r;
numTrialsPerLoc = D.numUnityTrialsPerLoc;
binSize         = 0.1;
x_range         = max(locResp(:)) - min(locResp(:));
x_min           = min(locResp(:)) - 0.2.*x_range; 
x_max           = max(locResp(:)) + 0.2.*x_range;
x_ticks         = s_V;
y_max           = 0.25; 
ybd             = [-y_max, y_max]; 
y_ticks         = linspace(-y_max,y_max,5);
modality        = {'A','V'};
modality_sign   = [1,-1];
modality_idx    = {'i','j'};
x               = x_min:binSize:x_max;
selected_plot   = [1]; %only plot the row that corresponds to s_A = 1
p_MAP           = cell(lenC, lenP);

for ii = 1:lenC
    for jj = 1:lenP
        %get predicted p(r|s_A, s_V)
        p_mAmV_mat      = R2{ii,jj}.p_mAmV_given_sAsV;
        MAP_mat         = round(squeeze(R2{ii,jj}.MAP),1);
        locResp_mat     = squeeze(locResp(ii,jj,:,:,:,:));
        [p_MAP{ii,jj}, counts] = predict_prob_r(MAP_mat, p_mAmV_mat,locResp_mat,...
                            x, sigma_r, model);
                     
        %plot
        if bool_plot
            for i = 1:lenS
                figure
                for j = 1:lenS
                    subplot(1, lenS, j)
                    addBackground([x_min, x_max], ybd,[x_min, x_ticks, x_max], y_ticks);
                    for n = 1:lenM
                        %fill the area under the distributions
                        patch([x, fliplr(x)], [modality_sign(n).*squeeze(p_MAP{ii,jj}(i,j,n,:))',...
                            zeros(1,length(x))],cMAP{n},'FaceAlpha',0.2,...
                            'EdgeColor',cMAP{n},'EdgeAlpha',0.7,'lineWidth',2); hold on;
                        %plot the MAP distributions
                        h(n) = plot(x, modality_sign(n).*squeeze(counts(i,j,n,:))./...
                            numTrialsPerLoc,'Color', cMAP{n},'lineWidth',1.5); hold on;
                    end
                    plot(x_min:binSize:x_max, zeros(1,length(x_min:binSize:x_max)),...
                        'Color', [0,0,0],'lineWidth',2,'lineStyle','-'); hold off;box off;
                    if j == 1
                        ylabel('$p(r|s_{AV_l,A},s_{AV_l,V})$','interpreter','latex'); 
                        yticks(y_ticks);yticklabels({'0.25','0.125','0','0.125','0.25'}); 
                    else; yticks([]); end
                    xticks(x_ticks(j)); xlabel('Stimulus location (dva)');
                    if i == 1 && j==1
                        legend([h(1), h(2)], {'Auditory','Visual'}, 'Location',...
                            'northeast','FontSize',20); legend boxoff;
                    end
                    xlim([x_min, x_max]); ylim(ybd); set(gca,'FontSize',20);%25
                end
                %add title        
                sgtitle(['Condition: ',model.cond{ii}, 'ruent ; Phase: ',...
                    model.phase{jj}, '; p_{C=1} = ', num2str(round(pC1(ii,jj),2))]); 
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.6]);%0.4
                set(gcf,'PaperUnits','centimeters','PaperSize',[65 30]);%17

                if ismember(i,selected_plot) && bool_save
                    saveFig('locResp', M, G, model.cond{ii},model.phase{jj});
                end
            end
        end
    end
end

function [p_MAP, counts] = predict_prob_r(MAP_mat, p_mAmV_mat, resp, s, sigma_r, model)
global lenM lenS
%This function generates predicted p(r|s_A, s_V) for all different
%combinations of s_A and s_V
%p_mAmV_mat: the joint probability of m_A and m_V given s_A and s_V
%            size = lenS x lenS x model.numBins_A x model.numBins_V
%MAP_mat   : the MAP estimates given an auditory and a visual stimulus
%            size = lenS x lenS x lenM x model.numBins_A x model.numBins_V
[p_r_given_shat,p_MAP,counts] = deal(zeros(lenS, lenS, lenM,length(s)));
for i = 1:lenS %for each auditory location
    for j = 1:lenS %for each visual location
        %get the joint likelihood and the matrix for MAP estimates given 
        %a selected visual and auditory location
        p_mAmV_temp = squeeze(p_mAmV_mat(i,j,:,:));
        MAP_temp    = cell(1,length(model.modality)); 
        for k = 1:lenM
            MAP_temp{k}     = round(squeeze(MAP_mat(i,j,k,:,:)),1);
            locResp_s       = round(squeeze(resp(i,j,k,:)),1);
            counts(i,j,k,:) = histcounts(locResp_s, ...
                                [s, s(end)+diff(s(1:2))]-diff(s(1:2))/2); 
        end
        %keep adding probability densities
        for l = 1:model.numBins_A 
            for m = 1:model.numBins_V
                for n = 1:lenM
                    %n = 1: get p(r_A|shat_A(l), shat_V(m))
                    %n = 2: get p(r_V|shat_A(l), shat_V(m))
                    p_r_given_shat(i,j,n,:) = norm_dst(s,MAP_temp{n}(l,m),sigma_r,0);
                    p_MAP(i,j,n,:) = p_MAP(i,j,n,:) +  p_r_given_shat(i,j,n,:).*p_mAmV_temp(l,m);
                end
            end
        end
    end
end

function p = norm_dst(x,mu,sigma,t)
p = 1/sqrt(2*pi*sigma^2).*exp(-(x-mu).^2./(2*sigma^2)) + t;

function r_samples = parametric_simulation(s, p_MAP_dist, n_samp, n_runs)
%first make sure the probability distribution sums up to 1
p_MAP_dist        = p_MAP_dist./sum(p_MAP_dist);
%compute the cumulative sum of the normalized probabilities
p_MAP_dist_cumsum = cumsum(p_MAP_dist);
%draw n samples from a uniform distribution
samp              = rand(n_runs,n_samp);
%find which bin do those random samples fall into
idx               = discretize(samp, [0, p_MAP_dist_cumsum]);
%find the corresponding location
r_samples         = s(idx);


%==========================================================================
%                         HELPING FUNCTIONS
%==========================================================================
%Add background (similar to Seaborn)
function addBackground(x_bd, y_bd, x_ticks, y_ticks)
patch([x_bd(1), x_bd(2), x_bd(2), x_bd(1)], [y_bd(1), y_bd(1), y_bd(2), y_bd(2)],...
    [234,234,242]./255, 'EdgeAlpha',0); hold on
if length(xticks) >=3 
    for i = 2:length(x_ticks)-1
        plot([x_ticks(i), x_ticks(i)], [y_bd(1), y_bd(2)], 'Color',...
            [1,1,1],'lineWidth', 0.5); hold on;
    end
end
if length(yticks) >= 3
    for j = 2:length(y_ticks)-1
        plot([x_bd(1), x_bd(2)],[y_ticks(j), y_ticks(j)], 'Color',...
            [1,1,1],'lineWidth', 0.5); hold on;
    end
end

%HELPING FUNCTION: Save the figure given a customized file name
function saveFig(name, M, G, cond, phase)
global sI; global bool_inc; 
if nargin < 4
    optional = '';
elseif nargin < 5
    optional = ['_',cond];
else
    optional = ['_',cond,'_', phase];
end
saveas(gcf, ['ModelPredictions_',name, '_', M, '_',G,'_unityJdg',...
        num2str(bool_inc(1)),'_locResp',num2str(bool_inc(2)),...
        optional,'_',sI], 'pdf'); 

