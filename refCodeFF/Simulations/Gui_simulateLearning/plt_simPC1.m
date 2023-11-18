function plt_simPC1(subjN, subjI, idx_subj, sesNum, cond_selected, cMAP, bool_save)
%load simulations
global param_name grid_name Cond
C          = load('dict_saveP_update_pC1_diffSigmaAV.mat', 'cell_savedP');
param_sim  = C.cell_savedP{1}(idx_subj,:);
grid_sim   = C.cell_savedP{end-1}(idx_subj,:);
pC1_sim    = C.cell_savedP{end}{idx_subj,1}{1}; 
%{1} is simulated pC1, {2} is simulated p(reporting C1)
%size: 2 (condition: 1st- congurent, 2nd-incongruent) x numBins x numBins
            
for i = 1:length(param_name)
    eval(['param.',param_name{i},'=',num2str(param_sim(i)),';']);
end

numBins = str2double(grid_sim{end});
for j = 1:4
    temp      = grid_sim{j};
    idx_comma = find(temp == ','); 
    temp_lb   = str2double(temp(2:(idx_comma-1)));
    temp_ub   = str2double(temp((idx_comma+1):(end-1)));
    eval([grid_name{j},' = linspace(temp_lb, temp_ub, numBins);']);
end

%find out which is the best-fitting model
addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                'ModelFitting/Fits']));
Mcmp        = load('Results_modelComparison_overall.mat', 'ModelComparison');
bestM       = Mcmp.ModelComparison{end}{idx_subj};
%e.g., 'SW-samePC1cond-strategyMAP_MA_strategyUnity_posteriorC1'
dash_idx    = find(bestM == '-');
M_name      = bestM((dash_idx(1)+1):(dash_idx(end)-1));
stratG      = bestM((dash_idx(end)+1):end);
seg_idx     = find(bestM == '_'); 
ds_locResp  = bestM((seg_idx(1)+1):(seg_idx(2)-1)); %'MS' or 'MA'
ds_unityJdg = bestM((seg_idx(3)+1):end); %measurements, MAP, posteriorC1

%load estimated pCommon
D           = load(['ModelFitting_updatePrior_',M_name,'_',stratG,'_sub',...
                num2str(subjN),'.mat'], 'ModelFitting');
%load the best-fit parameters and minNLL for each model
bestP    = D.ModelFitting{end}.P_f;
best_ct  = bestP(7);
if strcmp(M_name, 'fullModel')
    best_pC1_pre  = bestP([8,10]); best_pC1_post = bestP([9,11]);
elseif strcmp(M_name, 'samePC1pre')
    best_pC1_pre  = bestP([8,8]); best_pC1_post = bestP([9,10]);
else
    best_pC1_pre  = bestP([8,8]); best_pC1_post = bestP([9,9]);
end

% find the sim pC1 that's closest to estimated pC1
[best_sigma_AV_A,best_sigma_AV_V] = deal(NaN(1,2));
for i = 1:length(Cond) %2 conditions
    abs_diff = abs(squeeze(pC1_sim(i,:,:)) - best_pC1_post(i));
    [~, min_idx] = min(abs_diff(:));
    [min_idx_r, min_idx_c] = ind2sub([numBins, numBins], min_idx);
    best_sigma_AV_A(i) = eval([grid_name{1}(1:17),Cond{i},'(min_idx_r)']);
    best_sigma_AV_V(i) = eval([grid_name{2}(1:17),Cond{i},'(min_idx_c)']);
end

%load data from adaptation phase
nTT            = 160; 
empReportingC1 = NaN(1,length(Cond));
%get the data we need
for i = 1:length(Cond)
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
        'Adaptation v2/Data/',subjI]));
    if i == 1
        F = load(['Adaptation_',Cond{i}, '_sub', num2str(subjN),'_session',...
            num2str(sesNum(i)),'.mat'], 'Adaptation_data');
        D.AVpairs(i,1,:)  = F.Adaptation_data{4}.arrangedLocs_deg; %A loc
        D.AVpairs(i,3,:)  = F.Adaptation_data{3}.arrangedLocs_deg; %V loc
        D.AVpairs(i,2,:)  = zeros(1,nTT); D.AVpairs(i,4,:) = zeros(1,nTT);
        empReportingC1(i) = sum(F.Adaptation_data{end}.unity==1)/...
            sum(~isnan(F.Adaptation_data{end}.unity));
    else
        F = load(['Adaptation_',Cond{i}, '_sub', num2str(subjN),'_session',...
            num2str(sesNum(i)),'.mat'], ['Adaptation_',Cond{i},'_data']);
        D.AVpairs(i,1,:)  = F.Adaptation_incongruent_data{4}.arrangedLocs_deg; %A loc
        D.AVpairs(i,3,:)  = F.Adaptation_incongruent_data{3}.arrangedLocs_deg; %V loc
        D.AVpairs(i,2,:)  = -F.Adaptation_incongruent_data{3}.timing_relative./2; %A time relative
        D.AVpairs(i,4,:)  = F.Adaptation_incongruent_data{3}.timing_relative./2; %V time relative
        empReportingC1(i) = sum(F.Adaptation_incongruent_data{end}.unity==1)/...
            sum(~isnan(F.Adaptation_incongruent_data{end}.unity));
    end
end
D.totalTrials = size(D.AVpairs,3);
D.numSims     = 1e2;

%run simulations again
param.a_A            = bestP(1);
param.b_A            = bestP(2);
param.sigmaP_spatial = 100;
param.muP_spatial    = 0;
numSims              = 100;
sim_pC1_post         = NaN(length(Cond), numSims);
sim_pC1_allT         = NaN(length(Cond), nTT + 1, numSims);
[sim_pC1_mean, sim_pC1_allT_CI_lb, sim_pC1_allT_CI_ub] = deal(NaN(length(Cond), nTT+1));

for i = 1:length(Cond) %condition
    param.pCommon    = best_pC1_pre(i);
    param.sigma_spatial_AV_A = best_sigma_AV_A(i);
    param.sigma_spatial_AV_V = best_sigma_AV_V(i);
    param.alpha = eval(['param.alpha_pC1_', Cond{i}]);
    for j = 1:numSims
        [~, sim_pC1_post(i,j), sim_pC1_allT(i,:,j)] = ...
            simulateLearningPhase_update_pC1(...
            param, D, i, ds_unityJdg, best_ct, ds_locResp);
    end
    
    for k = 1:(nTT+1)
        sim_pC1_t = sort(squeeze(sim_pC1_allT(i,k,:)));
        sim_pC1_mean(i,k) = mean(sim_pC1_t);
        sim_pC1_allT_CI_lb(i,k) = sim_pC1_t(floor(numSims*0.025));
        sim_pC1_allT_CI_ub(i,k) = sim_pC1_t(floor(numSims*0.975));
    end
end

%plot it
x_bds = [1,nTT+1]; x_ticks=[x_bds(1), nTT/4, nTT/2,3*nTT/4, nTT];
y_bds = [min(sim_pC1_allT_CI_lb(cond_selected,:))-0.1,max(sim_pC1_allT_CI_lb(cond_selected,:))+0.1]; 
y_ticks = round(sort([best_pC1_pre(cond_selected),best_pC1_post(cond_selected)]),2); 

figure        
addBackground(x_bds, y_bds, x_ticks, [y_bds(1),y_ticks, y_bds(end)]);
plot(1:(nTT+1), zeros(1,nTT+1), '-', 'LineWidth',3,'Color',[0.5,0.5,0.5]); hold on
h2=plot(1:(nTT+1),sim_pC1_mean(cond_selected,:),'LineWidth',3,'Color',[0.4,0.4,0.4]); hold on;
patch([1:(nTT+1), (nTT+1):-1:1], [sim_pC1_allT_CI_lb(cond_selected,:),...
    fliplr(sim_pC1_allT_CI_ub(cond_selected,:))],cMAP,...
    'FaceAlpha',0.5,'EdgeColor','none','EdgeAlpha',1); hold on;
h1=scatter(1,best_pC1_pre(cond_selected),700,'d','MarkerEdgeColor','k',...
    'MarkerFaceColor', cMAP, 'MarkerFaceAlpha',0.5,'lineWidth',1);hold on;
scatter((nTT+1),best_pC1_post(cond_selected),700,'d','MarkerEdgeColor','k',...
    'MarkerFaceColor', cMAP,'MarkerFaceAlpha',0.5,'lineWidth',1);hold on;
text(2,y_bds(end)-0.02, subjI,'fontSize',20); hold off
box off; legend([h1,h2], {'Model estimates','Mean simulations'},'Location','southwest'); legend boxoff
xlim(x_bds); xticks(x_ticks); ylim(y_bds);yticks(y_ticks);
xlabel('Trial number (t)'); ylabel('p_{C=1}'); set(gca,'FontSize',25); 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.5]);
set(gcf,'PaperUnits','centimeters','PaperSize',[20 20]);
if bool_save; saveas(gcf, ['simPC1_with_modelFits_', subjI], 'pdf'); end


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
