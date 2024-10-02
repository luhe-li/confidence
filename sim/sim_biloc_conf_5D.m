function [bi_loc, bi_conf, conf_var, norm_var, est_var, criteria] = sim_biloc_conf_5D(...
    aA, bA, sigA, sigV1, sigV2, muP, sigP, sigConf, pCommon, fixP)

sim_d = fixP.sim_d;
sigMotor = fixP.sigMotor;
bi_nrep = fixP.bi_nrep;
n_sA = numel(fixP.sA);
sigVs = [sigV1, sigV2];

[shat, conf_var, norm_var] = deal(NaN(n_sA, n_sA, 2, 2, bi_nrep));
for cue = 1:numel(sigVs)

    sigV = sigVs(cue);

    %simulate measurements, which are drawn from Gaussian distributions
    % stochasticity starts here
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

    %compute the final location estimates if we assume model averaging.
    %Based on this strategy, the final location estimate is the sum of the
    %two intermediate location estimates, weighted by the corresponding
    %causal structure.
    %Eq. 4 in Wozny et al., 2010
    shat(:,:,1,cue,:) = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
    shat(:,:,2,cue,:) = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;

    % add motor noise
    bi_loc = randn(size(shat)).*sigMotor + shat;

    % calculate inverse of posterior variance as the proxy for confidence 
    % given each model
    switch sim_d
        case 1
            conf_var(:,:,1,cue,:)= repmat(1/JA + 1/JP, [n_sA, n_sA, 1, 1, bi_nrep]);
            conf_var(:,:,2,cue,:)= repmat(1/JV + 1/JP, [n_sA, n_sA, 1, 1, bi_nrep]);
        case 2
            conf_var(:,:,1:2,cue,:)= repmat(1/JV + 1/JA + 1/JP, [n_sA, n_sA, 1, 1, bi_nrep]);
            conf_var(:,:,1,cue,post_C1<0.5)  = 1/JA + 1/JP;
            conf_var(:,:,1,cue,post_C1<0.5)  = 1/JV + 1/JP;
        case 3
            conf_var(:,:,1,cue,:)= 1./(post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JA + 1/JP) + post_C1.* post_C2 .* (sHat_A_C2 - sHat_C1).^2);
            conf_var(:,:,2,cue,:)= 1./(post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JV + 1/JP) + post_C1.* post_C2 .* (sHat_V_C2 - sHat_C1).^2);
    end
end

% use a monotonic function to transfer confidencce variable to the unit of
% pixel
radius = alpha .* conf_var;




% 
points = expected_gain(bi_loc, radius, fixP);

%% utility function

function points = expected_gain(loc, conf, fixP)
% simulate probability density function assuming Gaussian
% loc   : localization (pixel)
% conf  : variance

% experimental parameters
x_axis = fixP.x_axis;
max_point = fixP.maxPoint;
min_point = fixP.minPoint;
elbow = fixP.elbow;
screenX = fixP.screenX;

pdf = normpdf(x_axis, loc, conf);


confRangeMax = min([loc - 1, screenX - loc],[], 5);
optRadius = NaN(length(loc),1);

for i = 1:length(estX)
    confRange = 0 : confRangeMax(i);
    costFun = max( maxPoint - confRange ./ elbow .* (maxPoint - minPoint) , minPoint);

    erPDF = myPDF(i,(estX(i)+1) : (estX(i) + confRangeMax(i))); % error pdf
    erCDF = [myPDF(i,estX(i)) , cumsum(erPDF) .* 2 + myPDF(i,estX(i))]; % error cdf, i.e. erf

    eGain = costFun .* erCDF;
    [~,optRadius(i)] = max(eGain);
    expGain(i).eGain = eGain;
end


end

