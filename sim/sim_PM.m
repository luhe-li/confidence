function [bi_loc, bi_conf, func] = sim_PM(aA, bA, sigA, sigV1, sigV2, muP, sigP, sigConf, pCommon, fixP)
% observer selects the intermediate posterior

sigMotor = fixP.sigMotor;
bi_nrep = fixP.bi_nrep;
n_sA = fixP.n_sA;
sigVs = [sigV1, sigV2];

[shat] = deal(NaN(n_sA, n_sA, 2, 2, bi_nrep));
post = NaN(n_sA, n_sA, 2, 2, bi_nrep, numel(fixP.center_axis));
for rel = 1:numel(sigVs)

    sigV = sigVs(rel);

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

    %compute the final location estimates if we assume probability matching.
    %Based on this strategy, the final location estimate is the integrated
    %one with a probability of post_C1, and is the segregated one with a
    %probability of post_C2.
    %Eq. 6 in Wozny et al., 2010
    random_sample = rand(size(post_C1));
    slc_indices = (post_C1 <= random_sample);
    sd_A(slc_indices) = sqrt(1/(1/JA+1/JP));
    sd_V(slc_indices) = sqrt(1/(1/JV+1/JP));
    shat_A(slc_indices) = sHat_A_C2(slc_indices);
    shat_V(slc_indices) = sHat_V_C2(slc_indices);

    % combine auditory and visual response
    shat(:,:,1,rel,:) = shat_A;
    shat(:,:,2,rel,:) = shat_V;

    for xx = 1:numel(fixP.center_axis)
        post(:,:,1,rel,:,xx) = normpdf(fixP.center_axis(xx), shat_A, sd_A);
        post(:,:,2,rel,:,xx) = normpdf(fixP.center_axis(xx), shat_V, sd_V);
    end

end

% %%%% debug
% 
% mu1 = 3.42;
% mu2 = -18;
% sd1 = 1;
% sd2 = 9.48;
% post1 = normpdf(fixP.center_axis, mu1, sd1);
% post2 = normpdf(fixP.center_axis, mu2, sd2);
% 
% [test_radius] = eGain_MAP([post1;post2], [mu1; mu2], fixP.maxScore, fixP.minScore, fixP.elbow, fixP.center_axis);
% 
% %%%%

% optimal radius given posterior and estimate
post_2d = reshape(post, [prod([n_sA, n_sA, 2, 2, bi_nrep]), numel(fixP.center_axis)]);
shat_1d = reshape(shat, [prod([n_sA, n_sA, 2, 2, bi_nrep]), 1]);
[opt_radius, ~, func]= eGain_MAP(post_2d, shat_1d, fixP.maxScore, fixP.minScore, fixP.elbow, fixP.center_axis);

% motor noise to location estimation
bi_loc = randn(size(shat)).*sigMotor + shat;

% adjustment noise to confidence radius
bi_conf = randn(size(opt_radius)).*sigConf + opt_radius;
bi_conf = reshape(bi_conf, [n_sA, n_sA, 2, 2, bi_nrep]);

check_plot=0;
if check_plot

    dims = [n_sA, n_sA, 2, 2, bi_nrep];

    % check posterior of all location combinations
    figure; hold on
    for aa = 1:4
        for vv = 1:4
            lin_index = sub2ind(dims, aa, vv, 1, 1, 1);
            plot(post_2d(lin_index, :))
        end
    end

    % check posterior of one location, across multiple trials
    figure; hold on
    for tt = 1:bi_nrep
        lin_index = sub2ind(dims, 2, 3, 1, 1, tt);
        plot(post_2d(lin_index, :))
    end

    % check cost function of one specific location combination
    figure; hold on
    lin_index = sub2ind(dims, 1, 3, 1, 1, 1); % change index here
    plot(func(lin_index).costFun);
    plot(func(lin_index).erCDF);
    plot(func(lin_index).gainFun);
    xlabel('Confidence unit')

end

end

%% check cost functions/posteriors
