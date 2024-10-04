function out = sim_MS(aA, bA,  sigV, sigA, sigP, sigC, pCommon, model)

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

% initiate all responses from intermediate posterior of a common cause
[sd_A, sd_V] = deal(repmat(sqrt(1/(1/JA+1/JV+1/JP)), [size(post_C1)]));
[shat_A, shat_V] = deal(sHat_C1);

% select the intermediate posterior of separate causes only if post_c1<=0.5
slc_indices = (post_C1 <= 0.5);
sd_A(slc_indices) = sqrt(1/(1/JA+1/JP));
sd_V(slc_indices) = sqrt(1/(1/JV+1/JP));
shat_A(slc_indices) = sHat_A_C2(slc_indices);
shat_V(slc_indices) = sHat_V_C2(slc_indices);

% combine auditory and visual response
shat(:,:,1,:) = shat_A;
shat(:,:,2,:) = shat_V;

% simulate posterior pdf for each trial using center coordinate
for xx = 1:numel(model.center_axis)
    post(:,:,1,:,xx) = normpdf(model.center_axis(xx), shat_A, sd_A);
    post(:,:,2,:,xx) = normpdf(model.center_axis(xx), shat_V, sd_V);
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