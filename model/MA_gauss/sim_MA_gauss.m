function out = sim_MA_gauss(aA, bA,  sigV, sigA, sigP, sigC, pCommon, model)

bi_nrep = model.bi_nrep;
n_sA = model.n_sA;
muP = model.muP;
sigMotor = model.sigMotor;

% simulate measurements, which are drawn from Gaussian distributions
% stimulus location combination x num of trial
mA2d    = randn(size(model.bi_sA',1), bi_nrep).*sigA + repmat((model.bi_sA' * aA + bA),1, bi_nrep);
mV2d    = randn(size(model.bi_sV',1), bi_nrep).*sigV + repmat(model.bi_sV', 1, bi_nrep);

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
shat(:,:,1,:) = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
shat(:,:,2,:) = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;

% approximate the posterior pdf by using a gaussian with the empirical mean
% and sd of this gaussian for each trial using center coordinate, where mu
% is , var is = w \sigma_1^2 + (1 - w) \sigma_2^2 + w (1 - w) (\mu_1 - \mu_2)^2

for xx = 1:numel(model.center_axis)
    post(:,:,1,:,xx) = normpdf(model.center_axis(xx), squeeze(shat(:,:,1,:)), ...
        post_C1/(1/JA+1/JV+1/JP) + post_C2/(1/JA+1/JP)+ post_C1 * post_C2 * (sHat_C1 - sHat_A_C2).^2);
    post(:,:,2,:,xx) = normpdf(model.center_axis(xx), squeeze(shat(:,:,2,:)), ...
        post_C1/(1/JA+1/JV+1/JP) + post_C2/(1/JV+1/JP)+ post_C1 * post_C2 * (sHat_C1 - sHat_V_C2).^2);
end

% optimal radius given posterior and estimate
post_2d = reshape(post, [prod([n_sA, n_sA, 2, bi_nrep]), numel(model.center_axis)]);
shat_1d = reshape(shat, [prod([n_sA, n_sA, 2, bi_nrep]), 1]);
[opt_radius, opt_gain, ~]= eGain_MAP(post_2d, shat_1d, model.maxScore, model.minScore, model.elbow, model.center_axis);

% motor noise to location estimation
out.bi_loc = randn(size(shat)).*sigMotor + shat;

% adjustment noise to confidence radius
bi_conf = randn(size(opt_radius)).*sigC + opt_radius;
out.bi_conf = reshape(bi_conf, [n_sA, n_sA, 2, bi_nrep]);

out.opt_gain = reshape(opt_gain,[n_sA, n_sA, 2, bi_nrep]);

end