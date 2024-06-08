function [loc, conf] = simUni(...
        aA, bA, sigA_uni, sigV1_uni, sigV2_uni, muP, sigP, sigC, c1, c2, c3, lapse, fixP)

% simulate for unimodal session

% loc
nrep = fixP.uni_nrep;
mA    = randn(1, nrep).*sigA_uni + (fixP.uni_sA' * aA + bA);
mV1    = randn(1, nrep).*sigV1_uni + fixP.uni_sV';
mV2    = randn(1, nrep).*sigV2_uni + fixP.uni_sV';

shatA = (mA*sigP^2 + muP*sigA_uni^2)./(sigA_uni^2 + sigP^2);
shatV1 = (mV1*sigP^2 + muP*sigV1_uni^2)./(sigV1_uni^2 + sigP^2);
shatV2 = (mV2*sigP^2 + muP*sigV2_uni^2)./(sigV2_uni^2 + sigP^2);

loc(1,:,:) = randn(size(shatA)).*fixP.sigMotor + shatA;
loc(2,:,:) = randn(size(shatV1)).*fixP.sigMotor + shatV1;
loc(3,:,:) = randn(size(shatV2)).*fixP.sigMotor + shatV2;

% conf
var(1,:,:) = repmat(1/(1/sigP^2 + 1/sigA_uni^2), [1, size(loc, 2, 3)]);
var(2,:,:) = repmat(1/(1/sigP^2 + 1/sigV1_uni^2), [1, size(loc, 2, 3)]);
var(3,:,:) = repmat(1/(1/sigP^2 + 1/sigV2_uni^2), [1, size(loc, 2, 3)]);
m = var;
v = sigC^2;
mu = log((m.^2)./sqrt(v+m.^2));
sigma = sqrt(log(v./(m.^2)+1));
est_var = lognrnd(mu, sigma);

% normalize
temp = reshape(est_var, 3, numel(fixP.uni_sA') * nrep);
mean_est_var = mean(temp, 2);
std_est_var = std(temp, [], 2);
norm_temp = (temp - mean_est_var)./ std_est_var;
norm_est_var = reshape(norm_temp, [3, numel(fixP.uni_sA'), nrep]);

% compare variance to criterion s.t. the lowest variance leads to highest
% confidence
conf = NaN(size(norm_est_var));
conf(norm_est_var < c1) = 4; 
conf(norm_est_var >= c1 & norm_est_var < c2) = 3;
conf(norm_est_var >= c2 & norm_est_var < c3) = 2;
conf(norm_est_var >= c3) = 1;

% add lapse to confidence report
lapse_trial = rand(size(conf))<lapse;

% Assign random confidence values between 1 and 4 for lapsed trials
conf(lapse_trial) = randi([1, 4], size(conf(lapse_trial)));

end