function [bi_loc, bi_conf, variance, norm_var, est_var] = simAllModels(...
    aA, bA, sigA, sigV, muP, sigP, pCommon, sigC, c1, c2, c3, lapse, fixP)

%% generative model of bimodal session

model_ind = fixP.model_ind;
sigMotor = fixP.sigMotor;
bi_nrep = fixP.bi_nrep;

%simulate measurements, which are drawn from Gaussian distributions
% stochasticity starts here
mA    = randn(1, bi_nrep).*sigA + (fixP.bi_sA * aA + bA);
mV    = randn(1, bi_nrep).*sigV + fixP.bi_sV;

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
shat(1,:) = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
shat(2,:) = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;
s
% add motor noise
bi_loc = randn(size(shat)).*sigMotor + shat;

% calculate posterior variance as the proxy for uncertainty (i.e., inverse
% of confidence) given each model
switch model_ind
    case 1
        variance(1,:)= repmat(1/(1/JA + 1/JP), [1, bi_nrep]);
        variance(2,:)= repmat(1/(1/JV + 1/JP), [1, bi_nrep]);
    case 2
        variance(1,:)= post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JA + 1/JP);
        variance(2,:)= post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JV + 1/JP);
%         % CONSIDER PREVIOUS M2 AS M3-MODEL SELECTION VER
%         variance(1:2,:)= repmat(1/(1/JV + 1/JA + 1/JP), [2, num_rep]);
%         variance(1,post_C1<0.5)  = 1/(1/JA + 1/JP);
%         variance(2,post_C1<0.5)  = 1/(1/JV + 1/JP);
    case 3
        variance(1,:)= post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JA + 1/JP) + post_C1.* post_C2 .* (sHat_A_C2 - sHat_C1).^2;
        variance(2,:)= post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JV + 1/JP) + post_C1.* post_C2 .* (sHat_V_C2 - sHat_C1).^2;
end

% make noisy measurements of variance/uncertainty for each modality
m = variance;
v = sigC^2;
mu = log((m.^2)./sqrt(v+m.^2));
sigma = sqrt(log(v./(m.^2)+1));
est_var = lognrnd(mu, sigma);

% normalize
min_est_var = min(est_var,2);
l_idx = nonzero(min_est_var - min(min_est_var));

s_idx = min(min(est_var,2));
est_var(s_idx,:)



mean_est_var = mean(est_var, 2);
std_est_var = std(est_var, [], 2);
norm_est_var = (est_var - mean_est_var)./ std_est_var;

% compare variance to criterion s.t. the lowest variance leads to highest
% confidence
bi_conf = NaN(size(norm_est_var));
bi_conf(norm_est_var < c1) = 4; 
bi_conf(norm_est_var >= c1 & norm_est_var < c2) = 3;
bi_conf(norm_est_var >= c2 & norm_est_var < c3) = 2;
bi_conf(norm_est_var >= c3) = 1;

% add lapse to confidence report
lapse_trial = rand(size(bi_loc))<lapse;

% Assign random confidence values between 1 and 4 for lapsed trials
bi_conf(lapse_trial) = randi([1, 4], size(bi_conf(lapse_trial)));

end