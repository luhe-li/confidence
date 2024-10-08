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
for xx = 1:numel(model.center_axis)
    post(:,1,:,xx) = normpdf(model.center_axis(xx), shatA, sqrt(1/(JA + JP)));
    post(:,2,:,xx) = normpdf(model.center_axis(xx), shatV, sqrt(1/(JV + JP)));
end
post_2d = reshape(post, [prod([n_sA, 2, n_rep]), numel(model.center_axis)]);
shat_1d = reshape(shat, [prod(n_sA, 2, n_rep), 1]);
[opt_radius, opt_gain, ~]= eGain_MAP(post_2d, shat_1d, model.maxScore, model.minScore, model.elbow, model.center_axis);

% location x post-cue/response x repetition
out.uni_loc = randn(size(shat)).*model.sigma_motor + shat;

% adjustment noise to confidence radius
uni_conf = randn(size(opt_radius)).*sigC + opt_radius;
out.uni_conf = reshape(uni_conf, [n_sA, 2, n_rep]);
out.opt_gain = reshape(opt_gain, [n_sA, 2, n_rep]);

end