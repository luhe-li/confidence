function out = sim_uni(aA, bA, sigma_A, sigma_V, sigma_P, mu_P, sigC, model)

% simulate measuremments, uni locations x num of trial
n_rep = model.uni_nrep;
n_sA = numel(model.uni_sA);
JP = 1/sigma_P^2;
JA = 1/sigma_A^2;
JV = 1/sigma_V^2;

mA    = randn(1, n_rep).*sigma_A + (model.uni_sA' * aA + bA);
mV    = randn(1, n_rep).*sigma_V + model.uni_sV';

shatA = (mA.*JA + mu_P.*JP)./(JA + JP);
shatV = (mV.*JV + mu_P.*JP)./(JV + JP);
shat(:,1,:) = shatA;
shat(:,2,:) = shatV;

% simulate posterior pdf for each trial using center coordinate
% only use 1 trial for each location because optimal confidence radius 
% only depends on the variance of posterior and elbow, is independent 
% of stimulus locations, and is the same across trials
% dimension of post_2d: 2(example trial for postcue each) x posterior axis
post_2d(1,:) = normpdf(model.center_axis, shatA(1,1), sqrt(1/(JA + JP)));
post_2d(2,:) = normpdf(model.center_axis, shatV(1,1), sqrt(1/(JV + JP)));
shat_1d = [shatA(1,1); shatV(1,1)];

[opt_radius, opt_gain, ~]= eGain_MAP(post_2d, shat_1d, model.maxScore, model.minScore, model.elbow, model.center_axis);

% location x post-cue/response x repetition
out.uni_loc = randn(size(shat)).*model.sigma_motor + shat;

% replicate radius for all locations and trials
expanded_radius = repmat(opt_radius', [n_sA, 1]);
uni_conf = repmat(expanded_radius, [1,1,n_rep]);
expanded_gain = repmat(opt_gain',[n_sA, 1]);
out.opt_gain = repmat(expanded_gain, [1,1,n_rep]);

% adjustment noise to confidence radius
out.uni_conf = randn(size(uni_conf)).*sigC + uni_conf;

end