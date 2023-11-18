%% subjects info
clear all; close all; clc
%subject info
subjNs  = [   3,   4,   5,   6,   8,   9,  11,  12,  13,  18,  19,  20,...
             21,  22,  23,  24,  25]; %15,16,17,
subjIs  = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','ZL','RE','MM',...
           'LL','CS','JH','MD','HHL'}; %outliers: 'AD','SM','SX'
lenS    = length(subjNs);

%%
[Prct_Outliers, sigma_r] = deal(NaN(1, lenS));
for i = 1:lenS
    addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
                            'Pointing task/Data/', subjIs{i}]));
    D  = load(['PointingTest_sub', num2str(subjNs(i)), '.mat'], 'PointingTest_data');
    data_pointing      = D.PointingTest_data{end}(1:2,:);
    cursorLoc          = data_pointing(1,:);
    cursor_locResp     = data_pointing(2,:);
    %remove outliers
    locError           = cursor_locResp - cursorLoc;
    SD                 = std(locError);
    bool_nonoutliers   = (locError./SD < 3 & locError./SD > -3);
    Prct_Outliers(i)   = sum(~bool_nonoutliers)/length(bool_nonoutliers);
    cursor_locResp_rm  = locError(bool_nonoutliers);
    
    %We use MLE to calculate sigma_r, which is sqrt(SSE/N)
    %We can't just use function std.m because the denominator is N-1
    sigma_r(i)         = sqrt(sum((cursor_locResp_rm - mean(cursor_locResp_rm)).^2)/...
                            length(cursor_locResp_rm));
end

%% print out the mean
fprintf('Mean sigma_r = %.4f, min = %.4f, max = %.4f. \n', mean(sigma_r),...
    min(sigma_r), max(sigma_r));

fprintf('Percentage of ourliers = %.4f, min = %.4f, max = %.4f. \n', mean(Prct_Outliers),...
    min(Prct_Outliers), max(Prct_Outliers));

%% save to a file
sigma_r = {subjNs, subjIs, sigma_r, [mean(sigma_r),...
    min(sigma_r), max(sigma_r)], [mean(Prct_Outliers),...
    min(Prct_Outliers), max(Prct_Outliers)]};
save('Results_sigma_r.mat','sigma_r');

