clear;
close all;

numA = 30;
sA = linspace(-30, 30, numA);
sV = sA;
sAV = combvec(sA, sV);
num_s = size(sAV,2);
num_trial = 100;

aA        = 1;
bA        = 0;

sigma_A  = 10;
sigma_V = 2.5;

mu_P = 0;
sigma_P = 33;

p_common = 0.8;

ds        = {'model averaging', 'model selection', 'probability matching'};
nDS       = length(ds);

bimResp = NaN(num_s, length(ds), 2, num_trial);

for i = 1:num_s

    % stimulus combination
    % 3 decision strategies (1. averaging, 2. selection, 3. matching)
    % 2 modalities (1. auditory, 2. visual)
    % num_trial
    bimResp(i,:,:,:) = simResp_CI_solutions(p_common, num_trial, sAV(1,i),...
        sAV(2,i), aA, bA, sigma_A, sigma_V, mu_P, sigma_P);

end

%% plot

d = 1;
m = 1;

mean_a_resp = mean(squeeze(bimResp(:, d, m, :)),2);
aResp = reshape(mean_a_resp, numA, numA);

figure; hold on
imagesc(aResp)
colorbar

xlabel('sV')
ylabel('sA')
% xticks(sV)
% xticklabels(sV)
% yticklabels(sA)
