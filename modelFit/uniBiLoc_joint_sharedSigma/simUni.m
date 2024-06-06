function [loc, conf] = simUni(...
        aA, bA, sigA_uni, sigV1_uni, sigV2_uni, muP, sigP, sigCA, sigCV, c1, c2, c3, lapse, fixP)

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
norm_var(1,:,:) = repmat(sigA_uni * sigP^2 / (sigP^2 + sigA_uni^2), [1, size(loc, 2, 3)]);
norm_var(2,:,:) = repmat(sigV1_uni * sigP^2 / (sigP^2 + sigV1_uni^2), [1, size(loc, 2, 3)]);
norm_var(3,:,:) = repmat(sigV2_uni * sigP^2 / (sigP^2 + sigV2_uni^2), [1, size(loc, 2, 3)]);
m = norm_var;
v(1,:,:) = repmat(sigCA,[1, size(norm_var, 2,3)]);
v(2:3,:,:) = repmat(sigCV,[2, size(norm_var, 2,3)]);
mu = log((m.^2)./sqrt(v+m.^2));
sigma = sqrt(log(v./(m.^2)+1));
est_var = lognrnd(mu, sigma);

% compare variance to criterion s.t. the lowest variance leads to highest
% confidence
conf = NaN(size(est_var));
conf(est_var < c1) = 4; 
conf(est_var >= c1 & est_var < c2) = 3;
conf(est_var >= c2 & est_var < c3) = 2;
conf(est_var >= c3) = 1;

% add lapse to confidence report
lapse_trial = rand(size(conf))<lapse;

% Assign random confidence values between 1 and 4 for lapsed trials
conf(lapse_trial) = randi([1, 4], size(conf(lapse_trial)));

end