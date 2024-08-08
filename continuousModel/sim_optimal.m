function [bi_loc, bi_conf] = sim_optimal(aA, bA, sigA, sigV1, sigV2, muP, sigP, sigConf, pCommon, fixP)

sigMotor = fixP.sigMotor;
bi_nrep = fixP.bi_nrep;
n_sA = fixP.n_sA;
sigVs = [sigV1, sigV2];

% [shat] = deal(NaN(n_sA, n_sA, 2, 2, bi_nrep));
for cue = 1:numel(sigVs)

    sigV = sigVs(cue);

    % simulate measurements, which are drawn from Gaussian distributions
    % stimulus location combination x num of trial
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

    % simulate posterior pdf for each trial using center coordinate
    for xx = 1:numel(fixP.center_axis)
        post(:,:,1,cue,:,xx) = post_C1.*normpdf(fixP.screen_cm(xx), sHat_C1, repmat(sqrt(const1), size(sHat_C1))) + post_C2.*normpdf(fixP.screen_cm(xx), sHat_A_C2, repmat(sqrt(constA), size(sHat_A_C2)));
        post(:,:,2,cue,:,xx) = post_C1.*normpdf(fixP.screen_cm(xx), sHat_C1, repmat(sqrt(const1), size(sHat_C1))) + post_C2.*normpdf(fixP.screen_cm(xx), sHat_V_C2, repmat(sqrt(constV), size(sHat_V_C2)));
    end

end

% for each possible estimate, given the posterior, calculate the optimal
% radius
post_2d = reshape(post, [prod([n_sA, n_sA, 2, 2, bi_nrep]), numel(fixP.screen_cm)]);

for xx = 1:numel(fixP.center_axis) 
    [radius(:,xx), gain(:,xx)] = eGain(post_2d, ones(size(post_2d,1),1).* xx, fixP.maxScore, fixP.minScore, fixP.elbow, fixP.center_axis;
end

% find max gain across all possible estimated pixels for each trial
[~, idx_opt] = max(gain, [], 2);
opt_est = fixP.px_axis(idx_opt) - 512; % convert back to -512:512 center coordinate
opt_radius = zeros(size(fixP.px_axis));
for tt = 1:numel(idx_opt)
    opt_radius(tt) = radius(tt,idx_opt(tt));    
end

% motor noise to location estimation
bi_loc = randn(size(opt_est)).*sigMotor + opt_est;
bi_loc = reshape(bi_loc, [n_sA, n_sA, 2, 2, bi_nrep]);

% adjustment noise to confidence radius
bi_conf = randn(size(opt_radius)).*sigConf + opt_radius;
bi_conf = reshape(bi_conf, [n_sA, n_sA, 2, 2, bi_nrep]);

end
