function [bi_loc, bi_conf, variance, norm_var, est_var] = simAllModels_5D_autoCriteria(...
    aA, bA, sigA, sigVs, muP, sigP, pCommon, sigC, lapse, fixP)

model_ind = fixP.model_ind;
sigMotor = fixP.sigMotor;
bi_nrep = fixP.bi_nrep;
n_sA = fixP.n_sA;

[shat, variance, norm_var] = deal(NaN(n_sA, n_sA, 2, 2, bi_nrep));

for cue = 1:numel(sigVs)

    sigV = sigVs(cue);

    %simulate measurements, which are drawn from Gaussian distributions
    % stochasticity starts here
    mA2d    = randn(size(fixP.bi_sA_comb',1), bi_nrep).*sigA + repmat((fixP.bi_sA_comb' * aA + bA),1, bi_nrep);
    mV2d    = randn(size(fixP.bi_sV_comb',1), bi_nrep).*sigV + repmat(fixP.bi_sV_comb', 1, bi_nrep);

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
        case 1 % suboptimal-select
            variance(:,:,1:2,cue,:) = repmat(1/(1/JV + 1/JA + 1/JP), [n_sA, n_sA, 2, 1, bi_nrep]);
            temp(post_C1<0.5) = 1/(1/JA + 1/JP);
            variance(:,:,1,cue,:) = reshape(temp, [n_sA, n_sA, bi_nrep]);
            temp(post_C1<0.5) = 1/(1/JV + 1/JP);
            variance(:,:,2,cue,:)  = reshape(temp, [n_sA, n_sA, bi_nrep]);
%             variance(1,post_C1<0.5)  = 1/(1/JA + 1/JP);
%             variance(2,post_C1<0.5)  = 1/(1/JV + 1/JP);
        case 2 % suboptimal-weight
            variance(:,:,1,cue,:)= post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JA + 1/JP);
            variance(:,:,2,cue,:)= post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JV + 1/JP);

        case 3 % optimal
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

num_c = 20;
autoc = linspace(min(est_var,[],'all'), max(est_var, [], 'all'), num_c);

bi_conf = NaN(size(est_var));
for ic = 2:num_c-1
    lc = autoc(ic-1);
    uc = autoc(ic);
    bi_conf(est_var>=lc & est_var<uc) = num_c - ic +1 ;
end

bi_conf(est_var<autoc(1)) = num_c;
bi_conf(est_var>=autoc(num_c)) = 1;

% add lapse to confidence report
lapse_trial = rand(size(bi_loc))<lapse;

% Assign random confidence values between 1 and 4 for lapsed trials
bi_conf(lapse_trial) = randi([1, num_c], size(bi_conf(lapse_trial)));

end