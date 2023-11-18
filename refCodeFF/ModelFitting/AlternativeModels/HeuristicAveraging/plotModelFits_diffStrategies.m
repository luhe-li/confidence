function [p, nLL] = plotModelFits_diffStrategies(modelTypes, stratG, bool_plot, bool_save)
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

%first load data, model info, best-fit p, minimum LL, sigma_r, s_A and s_V
global lenMT lenDS
[data, modelInfo_cell, p, nLL] = loadModelFits(modelTypes, stratG);
%get model predictions
for i = 1:lenMT
    for j = 1:lenDS
        modelInfo = modelInfo_cell{i,j};
        %calculate model predictions for each model
        %R1: mu_shat_A_uni, mu_shat_V_uni, sigma_shat_A_wN, sigma_shat_V_wN 
        %R2: pC1_given_sAsV, p_MAP, p_MAP_edges, counts
        %R3: s_A_prime, s_V_prime, s_A_hat, s_V_hat
        [~, R2, R3, pC1] = getModelPredictions(p{i,j}, modelTypes{i},...
            stratG{j}, data, modelInfo);
        if bool_plot(1) == 1 
            plotUnity(data.s_V, data, modelInfo, R2, pC1, modelTypes{i},...
                stratG{j}, bool_save(1));
        end
        if bool_plot(2) == 1 
            [loc_s,p_MAP] = plotlocResp(data.s_V, R2, R3, data, pC1, modelInfo,...
                modelTypes{i},stratG{j}, bool_save(2)); 
            VE_modelPred = predictVE(loc_s, p_MAP,[R3.shat_A;R3.shat_V]);
            if bool_plot(3) == 1
                plotVE(data.s_V, VE_modelPred,modelTypes{i},bool_save(3)); 
            end
        end
    end
end

%==========================================================================
%                        LOADING INFORMATION
%==========================================================================
function [data, MInfo, bestP, minNLL] = loadModelFits(M, stratG)
global sN sI lenMT lenDS
%add path and initialize a cell and a matrix
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'ModelFitting/Fits/', sI,'/Fits_alternativeModels/HeuristicAveraging']))
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
    'ModelFitting/Fits/', sI,'/Fits_locResp_conditioned_unityJdg_jointly']))
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


function [R1, R2, R3, pC1] = getModelPredictions(p, M, DS, data, modelInfo)
%calculate model prediction by passing the best-fit parameters to nLL funciton
if ~strcmp(DS, 'strategyMAP_HA_strategyUnity_MAP')
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
    end
else
    switch M
        case 'fullModel'
            [~, R1, R2] = nLL_overall_HA(p(1),p(2),p(3),p(4),p(5),p(6),p(7),...
                p(8),p(9),p(10),p(11),p(12),p(13),p(14),p(15),data,modelInfo);
            pC1 = round([p(8), p(9); p(10), p(11)],2);
        case 'samePC1pre'
            [~, R1, R2] = nLL_overall_samePC1pre_HA(p(1),p(2),p(3),p(4),p(5),...
                p(6),p(7),p(8),p(9),p(10),p(11),p(12),p(13),p(14),data,modelInfo);
            pC1 = round([p(8), p(9); p(8), p(10)],2);
        case 'samePC1cond'
            [~, R1, R2] = nLL_overall_samePC1cond_HA(p(1),p(2),p(3),p(4),p(5),...
                p(6),p(7),p(8),p(9),p(10),p(11),p(12),p(13),data,modelInfo);
            pC1 = round([p(8), p(8); p(9), p(9)],2);      
    end
end
%caculate s_A_prime, s_V_prime, s_V_hat
R3.s_A_prime = data.s_A.*p(1) + p(2);
R3.s_V_prime = data.s_V;
sigma_P      = 100;  mu_P = 0;
c_A          = (1/p(3)^2)/(1/p(3)^2+1/sigma_P^2); 
c_V          = (1/p(4)^2)/(1/p(4)^2+1/sigma_P^2);
f_A          = (mu_P/sigma_P^2)/(1/p(3)^2+1/sigma_P^2);
f_V          = (mu_P/sigma_P^2)/(1/p(4)^2+1/sigma_P^2);
R3.shat_V    = R3.s_V_prime.*c_V + f_V; R3.shat_A = R3.s_A_prime.*c_A + f_A;

function VE_modelPred = predictVE(s, p_MAP, shat)
global lenD lenS lenC lenP lenM
%size(p_MAP): lenC x lenP (length(s_A) x length(s_V) x lenM x length(s))
VE_modelPred = NaN(lenC, lenP, lenM, lenD);
for i = 1:lenC
    for j = 1:lenP
        p_MAP_mat = p_MAP{i,j};
        [VE_temp, inv_VE_temp] = deal(NaN(lenS,lenD));
        for k = 1:lenS %V loc 
            %by model predictions
            for l = 1:lenS %A loc
                %Auditory loc shifts
                p_MAP_norm = squeeze(p_MAP_mat(l,k,1,:))';
                p_MAP_norm = p_MAP_norm./sum(p_MAP_norm);
                VE_temp(k,ceil(lenD/2)+l-k) = shat(2,l) - p_MAP_norm*s';
                %Visual loc shifts
                p_MAP_norm = squeeze(p_MAP_mat(k,l,2,:))';
                p_MAP_norm = p_MAP_norm./sum(p_MAP_norm);
                inv_VE_temp(k,ceil(lenD/2)+l-k) = p_MAP_norm*s' - shat(1,l);
            end
        end
        VE_modelPred(i,j,1,:) = nanmean(VE_temp);
        VE_modelPred(i,j,2,:) = nanmean(inv_VE_temp);
    end
end

%==========================================================================
%                         PLOTTING FUNCTIONS
%==========================================================================
function plotVE(s_V, VE_modelPred, M, bool_save)
%shat_V and shat_A are the perceived stimulus location 
global sI sN lenM lenC lenP lenD; 
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                'Bimodal localization v2/Data/',sI]));
C              = load(['Summary_VE_sub', num2str(sN),'.mat'], 'Summary_VE');
meanVE         = C.Summary_VE{2};
sdLocR         =  C.Summary_VE{3};
numT_AV        = [20,40,60,80,60,40,20];
AV_combs       = combvec(s_V, s_V);
AV_offset      = sort(unique(AV_combs(2,:) - AV_combs(1,:)));
modality       = {'A','V'};
cond           = {'Congruent', 'Incongruent'};
y_lbl          = {'Auditory shifts (deg)', 'Visual shifts (deg)'};
lgd_cond_phase = {'Pre-learning', 'Post-learning'};
x_bds_VE       = [AV_offset(1)-3, AV_offset(end)+3];
x_ticks_VE     = [x_bds_VE(1), AV_offset, x_bds_VE(end)];
y_bds_VE       = [AV_offset(1)-2, AV_offset(end)+2];
y_ticks_VE     = AV_offset(find(y_bds_VE(1) < AV_offset,1):...
                           find(y_bds_VE(end) > AV_offset,1,'last'));
cb             = [65,105,225]./255; cr = [0.85, 0.23, 0.28];
cMap_VE        = {min(cb.*2,1),cb; [247,191,190]./255, cr};
jitter         = [-0.8,0.8];
lw             = 3; %lineWidth
fs_lbls        = 25;
fs_lgds        = 20;

for i = 1:lenM
    figure
    for j = 1:lenC
        subplot(1,lenC,j)
        addBackground(x_bds_VE, y_bds_VE, x_ticks_VE, [y_bds_VE(1),y_ticks_VE,y_bds_VE(end)])
        if i == 1
            plt3 = plot([AV_offset(1),AV_offset(end)],[AV_offset(1),AV_offset(end)],...
                'k--','Color',ones(1,3).*0.5,'lineWidth',lw); hold on;
        else
            plt3 = plot([AV_offset(1),AV_offset(end)],[AV_offset(end),AV_offset(1)],...
                'k--','Color',ones(1,3).*0.5,'lineWidth',lw); hold on;
        end
        plt4 = plot(AV_offset, zeros(1,lenD), 'k:','Color',...
            ones(1,3).*0.5,'lineWidth',lw); hold on;
        for k = 1:lenP
            for n = 1:lenD
                plt(k) = errorbar(AV_offset(n)+jitter(k), meanVE(j,k,i,n),...
                    sdLocR(j,k,i,n), '-o','MarkerSize',sqrt(numT_AV(n)).*3,...
                    'Color',cMap_VE{i,k},'MarkerFaceColor',cMap_VE{i,k},...
                    'MarkerEdgeColor',cMap_VE{i,k},'lineWidth',lw); hold on;
            end
            plot(AV_offset+jitter(k), squeeze(VE_modelPred(j,k,i,:)),'-',...
                'lineWidth',lw,'Color',cMap_VE{i,k}); hold on
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


function plotUnity(s_V, data, model, R2, pC1, M, G, bool_save)
global sI lenC lenP lenS lenD
calcSD    = @(p,n) sqrt(p.*(1-p)./n);
numT_AV   = [20,40,60,80,60,40,20];
AV_combs  = combvec(s_V, s_V);
AV_offset = sort(unique(AV_combs(2,:) - AV_combs(1,:)));
x_bds     = [AV_offset(1) - 2, AV_offset(end) + 2]; y_bds = [-0.05,1.05];
x_ticks   = AV_offset; y_ticks = 0:0.25:1;
colorMap  = [0.65, 0.65,0.65;0.1,0.1,0.1];%[135,206,235; 65,105,225]./255;
jitter    = [-0.8,0.8];
lw        = 3; %lineWidth
fs_lbls   = 25;
fs_lgds   = 20;
lgd_cond_phase = {'Pre-learning', 'Post-learning'};

figure
for i = 1:lenC
    subplot(1,2,i)
    addBackground(x_bds, y_bds, x_ticks, y_ticks)
    %initialize the proportion of reporting C=1
    [overall_pC1_pre, overall_pC1_post,overall_pC1_pre_model, overall_pC1_post_model] =...
        deal(NaN(lenS,lenD));
    for k = 1:lenS
        idx_lb = ceil(lenD/2)+1-k;
        idx_ub = lenD + 1 -k;
        %by data
        overall_pC1_pre(k,idx_lb:idx_ub) = squeeze(data.bimodal_unity_prob(i,1,:,k))';
        overall_pC1_post(k,idx_lb:idx_ub) = squeeze(data.bimodal_unity_prob(i,2,:,k))';
        %by model predictions
        overall_pC1_pre_model(k,idx_lb:idx_ub) = squeeze(R2{i,1}.pC1_given_sAsV(:,k))';
        overall_pC1_post_model(k,idx_lb:idx_ub) = squeeze(R2{i,2}.pC1_given_sAsV(:,k))';
    end
    %plot the model predictions
    for j = 1:lenP
        plot(AV_offset, eval(['nanmean(overall_pC1_', model.phase{j},...
            '_model)']),'Color',colorMap(j,:),'lineWidth', lw, 'lineStyle','-'); hold on;
    end
    
    mean_lenD_pre  = nanmean(overall_pC1_pre);
    mean_lenD_post = nanmean(overall_pC1_post);
    SD_pre         = calcSD(mean_lenD_pre, numT_AV); 
    SD_post        = calcSD(mean_lenD_post, numT_AV); 
    for m = 1:length(AV_offset)
        plt_1 = errorbar(AV_offset(m)+jitter(1), mean_lenD_pre(m),...
            SD_pre(m),'-o','MarkerSize',sqrt(numT_AV(m)).*3,'Color',...
            colorMap(1,:),'MarkerFaceColor',colorMap(1,:),'MarkerEdgeColor',...
            colorMap(1,:),'lineWidth',lw); hold on
        plt_2 = errorbar(AV_offset(m)+jitter(2),mean_lenD_post(m),...
            SD_post(m),'-o','MarkerSize',sqrt(numT_AV(m)).*3,'Color',...
            colorMap(2,:),'MarkerFaceColor',colorMap(2,:),'MarkerEdgeColor',...
            colorMap(2,:),'lineWidth',lw); hold on
    end
%     %plot the behavioral data
%     plt_1 = scatter(AV_offset+jitter(1), nanmean(overall_pC1_pre), markerS,...
%         colorMap(1,:),'o','MarkerFaceColor',colorMap(1,:), ...
%         'MarkerFaceAlpha',0.5); hold on
%     plt_2 = scatter(AV_offset+jitter(2), nanmean(overall_pC1_post), markerS,...
%         colorMap(2,:),'d','MarkerFaceColor',colorMap(2,:), ...
%         'MarkerFaceAlpha',0.5); hold on;
    %add legends
    if i == 1; text(x_bds(1) + 0.2, 1.02, sI,'FontSize',fs_lgds); end
    if i == 1; legend([plt_1 plt_2], lgd_cond_phase,...
            'Position',[0.3 0.15 0.05 0.1]); legend boxoff; end 
    xticks(x_ticks); xlim(x_bds); xlabel('Spatial discrepancy (V - A, deg)'); 
    ylabel(sprintf(['The probability of \nreporting a common cause'])); 
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


function [x,p_MAP] = plotlocResp(s_V, R2, R3, D, pC1, model, M, G, bool_save)
global lenC lenP lenM lenS
cb = [65,105,225]./255; cr = [0.85, 0.23, 0.28];
cMAP = {cb, cr};
locResp         = D.bimodal_locResp;
sigma_r         = D.sigma_r;
numTrialsPerLoc = D.numUnityTrialsPerLoc;
binSize         = 0.1;
x_min           = -25; 
x_max           = 25;
x_ticks         = s_V;
y_max           = 0.35; 
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
                    %plot the mean of biased estimates
                    plot([eval(['R3.s_', modality{n},'_prime(',...
                        modality_idx{n},')']),eval(['R3.s_',modality{n},...
                        '_prime(', modality_idx{n},')'])],...
                        sort([0, modality_sign(n).*y_max]),...
                        'Color',cMAP{n},'lineWidth',1.5,'lineStyle','--'); hold on;
                end
                plot(x_min:binSize:x_max, zeros(1,length(x_min:binSize:x_max)),...
                    'Color', [0,0,0],'lineWidth',2,'lineStyle','-'); hold off;box off;
                if j == 1
                    ylabel('p(r|s_{AV,A},s_{AV,V})'); 
                    yticks(y_ticks);yticklabels({'0.25','0.125','0','0.125','0.25'}); 
                else; yticks([]); end
                xticks(x_ticks(j)); xlabel('Stimulus location (dva)');
                if i == 1 && j==1
                    legend([h(1), h(2)], {'Auditory','Visual'}, 'Location',...
                        'northeast','FontSize',20); legend boxoff;
                end
                xlim([x_min, x_max]); ylim(ybd); set(gca,'FontSize',25);
            end
            %add title        
            sgtitle(['Condition: ',model.cond{ii}, 'ruent ; Phase: ',...
                model.phase{jj}, '; p_{C=1} = ', num2str(round(pC1(ii,jj),2))]); 
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.4]);
            set(gcf,'PaperUnits','centimeters','PaperSize',[65 17]);

            if ismember(i,selected_plot) && bool_save
                saveFig('locResp', M, G, model.cond{ii},model.phase{jj});
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
                    p_MAP(i,j,n,:) = p_MAP(i,j,n,:) + ...
                        p_r_given_shat(i,j,n,:).*p_mAmV_temp(l,m);
                end
            end
        end
    end
end

function p = norm_dst(x,mu,sigma,t)
p = 1/sqrt(2*pi*sigma^2).*exp(-(x-mu).^2./(2*sigma^2)) + t;


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

