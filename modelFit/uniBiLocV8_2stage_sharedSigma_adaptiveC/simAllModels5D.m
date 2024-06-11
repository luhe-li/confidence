function [bi_loc, bi_conf, variance, norm_var, est_var, criteria] = simAllModels5D(...
    aA, bA, sigA, sigVs, muP, sigP, pCommon, sigC, c1, c2, c3, lapse, fixP)

%% generative model of bimodal session

model_ind = fixP.model_ind;
sigMotor = fixP.sigMotor;
bi_nrep = fixP.bi_nrep;
n_sA = 4;

[shat, variance, norm_var] = deal(NaN(n_sA, n_sA, 2, 2, bi_nrep));

for cue = 1:numel(sigVs)

    sigV = sigVs(cue);
    %simulate measurements, which are drawn from Gaussian distributions
    % stochasticity starts here
    mA2d    = randn(size(fixP.bi_sA',1), bi_nrep).*sigA + repmat((fixP.bi_sA' * aA + bA),1, bi_nrep);
    mV2d    = randn(size(fixP.bi_sV',1), bi_nrep).*sigV + repmat(fixP.bi_sV', 1, bi_nrep);

    mA = reshape(mA2d, [n_sA, n_sA, bi_nrep]);
    mV = reshape(mV2d, [n_sA, n_sA, bi_nrep]);

    %compute constants (these will come in handy when writing the equations
    %for the likelihood/posterior of a common cause and separate causes.
    JA    = sigA^2;
    JV    = sigV^2;
    JP    = sigP^2;
    const1= JA*JV + JA*JP + JV*JP;
    constA= JA + JP;
    constV= JV + JP;

    %calculate the likelihood of a common cause and separate causes
    %Eq. 4 in Körding et al., 2007
    L_C1  = 1/(2*pi*sqrt(const1)).*exp(-0.5.*((mA - mV).^2 .* ...
        JP + (mV - muP).^2.*JA + (mA - muP).^2.* JV)./const1);
    %Eq. 6 in Körding et al., 2007
    L_C2  = 1/(2*pi*sqrt(constA*constV)).*exp(-0.5.*(mA - muP).^2./constA +...
        (mV - muP).^2./constV);

    %calculate posterior of a common cause and separate causes
    %Eq. 2 in Körding et al., 2007
    post_C1    = pCommon.*L_C1./(pCommon.*L_C1 + (1-pCommon).*L_C2);
    %posterior of separate causes
    post_C2    = 1 - post_C1;

    %compute the two intermediate location estimates
    %An integrated intermediate estimate is the sum of mA, mV and muP with
    %each weighted by their relative reliabilities
    %Eq. 12 in Körding et al., 2007
    sHat_C1    = (mA./JA + mV./JV + muP/JP)./(1/JV + 1/JA + 1/JP);

    %A segregated intermediate estimate is the sum of mA/mV and muP with
    %each weighted by their relative reliabilities
    %Eq. 11 in Körding et al., 2007
    sHat_A_C2  = (mA./JA + muP/JP)./(1/JA + 1/JP);
    sHat_V_C2  = (mV./JV + muP/JP)./(1/JV + 1/JP);

    %compute the final location estimates if we assume model averaging.
    %Based on this strategy, the final location estimate is the sum of the
    %two intermediate location estimates, weighted by the corresponding
    %causal structure.
    %Eq. 4 in Wozny et al., 2010
    shat(:,:,1,cue,:) = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
    shat(:,:,2,cue,:) = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;

    % add motor noise
    bi_loc = randn(size(shat)).*sigMotor + shat;

    % calculate posterior variance as the proxy for uncertainty (i.e., inverse
    % of confidence) given each model
    switch model_ind
        case 1
            variance(:,:,1,cue,:)= repmat(1/(1/JA + 1/JP), [n_sA, n_sA, 1, 1, bi_nrep]);
            variance(:,:,2,cue,:)= repmat(1/(1/JV + 1/JP), [n_sA, n_sA, 1, 1, bi_nrep]);
        case 2
            variance(:,:,1,cue,:)= post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JA + 1/JP);
            variance(:,:,2,cue,:)= post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JV + 1/JP);
            %         % CONSIDER PREVIOUS M2 AS M3-MODEL SELECTION VER
            %         variance(1:2,:)= repmat(1/(1/JV + 1/JA + 1/JP), [2, num_rep]);
            %         variance(1,post_C1<0.5)  = 1/(1/JA + 1/JP);
            %         variance(2,post_C1<0.5)  = 1/(1/JV + 1/JP);
        case 3
            variance(:,:,1,cue,:)= post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JA + 1/JP) + post_C1.* post_C2 .* (sHat_A_C2 - sHat_C1).^2;
            variance(:,:,2,cue,:)= post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JV + 1/JP) + post_C1.* post_C2 .* (sHat_V_C2 - sHat_C1).^2;
    end

    % normalize variance by the corresponding modality noise
    norm_var(:,:,1,cue,:) = variance(:,:,1,cue,:)./sigA;
    norm_var(:,:,2,cue,:) = variance(:,:,2,cue,:)./sigV;

end

% make noisy measurements of variance/uncertainty for each modality
m = norm_var;
v = sigC;
mu = log((m.^2)./sqrt(v+m.^2));
sigma = sqrt(log(v./(m.^2)+1));
est_var = lognrnd(mu, sigma);

% find criteria based on quantile
quantiles = [c1, c2, c3];
all_criteria = quantile(est_var, quantiles, 5);

% 2d version: seprate for modality only
criteria = squeeze(mean(mean(mean(all_criteria, 4),2),1));
bi_conf = NaN(size(est_var));
for cue = 1:2

    q1 = criteria(cue, 1);
    q2 = criteria(cue, 2);
    q3 = criteria(cue, 3);

    i_est_var = squeeze(est_var(:, :, cue, :, :));
    temp_conf = NaN(size(i_est_var));
    temp_conf(i_est_var < q1) = 4;
    temp_conf(i_est_var >= q1 & i_est_var < q2) = 3;
    temp_conf(i_est_var >= q2 & i_est_var < q3) = 2;
    temp_conf(i_est_var >= q3) = 1;

    bi_conf(:, :, cue, :, :) = temp_conf;

end

% % 3d version: seprate for modality and reliability -> overlapping curves
% criteria = squeeze(mean(mean(all_criteria, 2),1));
% bi_conf = NaN(size(est_var));
% for cue = 1:2
% 
%     for rel = 1:2
%     q1 = criteria(cue, rel, 1);
%     q2 = criteria(cue, rel, 2);
%     q3 = criteria(cue, rel, 3);
% 
%     i_est_var = squeeze(est_var(:, :, cue, rel, :));
%     temp_conf = NaN(size(i_est_var));
%     temp_conf(i_est_var < q1) = 4;
%     temp_conf(i_est_var >= q1 & i_est_var < q2) = 3;
%     temp_conf(i_est_var >= q2 & i_est_var < q3) = 2;
%     temp_conf(i_est_var >= q3) = 1;
% 
%     bi_conf(:, :, cue, rel, :) = temp_conf;
% 
% end

% % 5d version: seprate for modality, reliability, and location -> flat
% for aa = 1:4
%     for vv = 1:4
%         for cue = 1:2
%             for rel = 1:2
% 
%                 q1 = criteria(aa, vv, cue, rel, 1);
%                 q2 = criteria(aa, vv, cue, rel, 2);
%                 q3 = criteria(aa, vv, cue, rel, 3);
% 
%                 i_est_var = squeeze(est_var(aa, vv, cue, rel, :));
%                 temp_conf = NaN(size(i_est_var));
%                 temp_conf(i_est_var < q1) = 4;
%                 temp_conf(i_est_var >= q1 & i_est_var < q2) = 3;
%                 temp_conf(i_est_var >= q2 & i_est_var < q3) = 2;
%                 temp_conf(i_est_var >= q3) = 1;
% 
%                 bi_conf(aa, vv, cue, rel, :) = temp_conf;
% 
%             end
%         end
%     end
% end

% visualize criteria to see if they are resonable
for cue = 1:2
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    tiledlayout(4,4)
    for sA = 1:4
        for sV = 1:4
            nexttile
            hold on

            for rel = 1:2
                histogram(squeeze(est_var(sA, sV, cue, rel, :)),'NumBins',25)
            end

            for cc = 1:3
                xline(criteria(cue, cc))
            end
        end
    end
end

% add lapse to confidence report
lapse_trial = rand(size(bi_loc))<lapse;

% Assign random confidence values between 1 and 4 for lapsed trials
bi_conf(lapse_trial) = randi([1, 4], size(bi_conf(lapse_trial)));

end