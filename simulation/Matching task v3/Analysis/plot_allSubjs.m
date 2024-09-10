clear all; close all; clc
%load the data
subjN_dict = [   3,   4,   5,   6,   8,   9,  11,  12,  13,  18,  19,  20,...
                21,  22,  23,  24,  25];%15AD, 16SM, 17SX are identified as outliers
subjI_dict = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','ZL','RE','MM',...
              'LL','CS','JH','MD','HHL'};
lenS       = length(subjN_dict);
get95CI    = @(m,n) [m(ceil(n*0.025)), m(floor(n*0.975))];
[a_A, b_A, lb_a_A, ub_a_A, lb_b_A, ub_b_A] = deal(NaN(1,lenS));
pltInfo.numBins   = 12; 
pltInfo.bool_save = 0;

[data, binnedD] = deal(cell(1, lenS));
for i = 1:lenS
    subjN       = subjN_dict(i); 
    subjI       = subjI_dict{i}; pltInfo.subjI = subjI;
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                     'Matching task v3/Data/', subjI]));
    C           = load(['AV_alignment_sub', num2str(subjN), '_dataSummary.mat'],...
                    'AV_alignment_data');
    ExpInfo     = C.AV_alignment_data{1};
    data{i}     = C.AV_alignment_data{2};
    numBtst     = size(data{i}.PSE,1);
    a_A(i)      = data{i}.polyfit(1); b_A(i) = data{i}.polyfit(2);
    PSE_btst    = data{i}.PSE;
    ab_btst     = arrayfun(@(idx) polyfit([-12,-4,4,12], PSE_btst(idx,:),1)',...
                        1:numBtst, 'UniformOutput', false);
    ab_btst_mat = [ab_btst{:}]; 
    bds_a_A     = get95CI(sort(ab_btst_mat(1,:)), numBtst); 
    lb_a_A(i)   = bds_a_A(1); ub_a_A(i) = bds_a_A(end);
    bds_b_A     = get95CI(sort(ab_btst_mat(2,:)), numBtst);
    lb_b_A(i)   = bds_b_A(1); ub_b_A(i) = bds_b_A(end);
    
    D               = load(['A_aligns_V_sub', num2str(subjN), '.mat'],'A_aligns_V_data');
    Distance        = D.A_aligns_V_data{8};
    Distance(:,end) = []; %the last one was generated but never used in the experiment
    %specify plotting information
    bds             = [min(Distance(:)), max(Distance(:))]; numX = 1e3;
    pltInfo.x       = linspace(bds(1)-5,bds(end)+5, numX);
    %plot psychometric curve
    binnedD{i}      = plot_psychfunc(data{i}, ExpInfo, pltInfo);
    
    %plot PSE
    plot_PSE(data{i}, ExpInfo, pltInfo)
end

%% compute R squared
cnorm = @(t,p) normcdf(t,p(1),p(2)).*(1-p(3))+p(3)/2;
R2_adjusted = NaN(1, lenS);
for i = 1:lenS
    data_i    = data{i};
    binnedD_i = binnedD{i};
    [SSE_i, SST_i, numDataPoints_i] = deal(NaN(1,ExpInfo.testLocations));
    for j = 1:ExpInfo.testLocations
        binnedD_ij      = binnedD_i{j};
        %binned x values
        binnedD_ij_x    = binnedD_ij(1,:);
        %percent of reporting A to the right (binned)
        binnedD_ij_prct = binnedD_ij(3,:)./binnedD_ij(2,:);
        %prediction by the best-fitting psychometric curves
        prct_prediction = cnorm(binnedD_ij_x,[data_i.estimatedP(j*2:(j*2+1)),...
            data_i.estimatedP(1)]);
        %Sum squared residuals
        SSE_i(j) = sum((binnedD_ij_prct - prct_prediction).^2);
        %Sum squared total
        SST_i(j) = sum((binnedD_ij_prct - mean(binnedD_ij_prct)).^2);
        numDataPoints_i(j) = length(binnedD_ij_x);
    end
    SSE            = sum(SSE_i);
    SST            = sum(SST_i);
    R2             = 1 - SSE/SST;
    %number of total datapoints
    n              = sum(numDataPoints_i);
    %number of free parameters
    k              = length(data_i.estimatedP);
    %adjusted R squared
    R2_adjusted(i) = 1 - (1-R2)*(n-1)/(n-k-1);
end
fprintf('Psychometric functions: Mean adjusted R2 = %.3f, min = %.3f, max = %.3f.\n', ...
    mean(R2_adjusted), min(R2_adjusted), max(R2_adjusted));
fprintf('Linear regressions: Mean R2 = %.3f, min = %.3f, max = %.3f.\n',...
    mean(arrayfun(@(idx) data{idx}.R_square, 1:lenS)), ...
    min(arrayfun(@(idx) data{idx}.R_square, 1:lenS)),...
    max(arrayfun(@(idx) data{idx}.R_square, 1:lenS)));

%% to answer reviewer 2's question
%"I wonder whether the relation between the PSEs and the location of the 
%standard visual stimulus is indeed a linear one; for instance, looking a 
%Figure 2B (right panel) one may also consider such a relation as locally 
%linear with a slope of 1 (plus an offset) between the two central points, 
%and shallower towards the two extreme points."
R2_locallyLinear = NaN(1, lenS);
for i = 1:lenS
    %first get all the four PSEs
    PSEs    = data{i}.estimatedP(2:2:end);
    %pass each participant's PSEs into function sqErr_locallyLinearFunc to
    %compute squared error
    errFunc = @(p) sqErr_locallyLinearFunc(p(1),p(2),-7.5:5:7.5, PSEs);
    %use fminsearch to find the best slope and intercept
    estimatedP_locallyLinear = fminsearch(errFunc, [rand, rand]);
    %compute R squared
    [SSE_locallyLinear, y_pred] = sqErr_locallyLinearFunc(estimatedP_locallyLinear(1), ...
        estimatedP_locallyLinear(2), -7.5:5:7.5, PSEs);
    %compute R2
    SST_locallyLinear = sum((PSEs - mean(PSEs)).^2);
    R2_locallyLinear(i) = 1 - SSE_locallyLinear/SST_locallyLinear;
end
fprintf('Locally linear functions: Mean R2 = %.3f, min = %.3f, max = %.3f.\n',...
    mean(R2_locallyLinear), min(R2_locallyLinear), max(R2_locallyLinear));

%% plot
mean_a_A = mean(a_A); 
mean_b_A = mean(b_A); 
SE_a_A   = std(a_A)/sqrt(lenS); 
SE_b_A   = std(b_A)/sqrt(lenS);
disp(['Group mean a_A: ', num2str(mean_a_A), ', SEM: ', num2str(SE_a_A)]);
disp(['Group mean b_A: ', num2str(mean_b_A), ', SEM: ', num2str(SE_b_A)]);
bool_save = 1;
plot_regLines(a_A, b_A, lb_a_A, ub_a_A, lb_b_A, ub_b_A, mean_a_A,...
    mean_b_A, SE_a_A, SE_b_A, bool_save)

