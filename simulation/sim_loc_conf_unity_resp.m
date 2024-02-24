function [loc, conf, unity, pdf] = sim_loc_conf_unity_resp(pCommon, num_rep, sA, sV, aA, bA, sigA, sigV, muP, sigP, fixP)

%SIM_LOC_CONF_UNITY_RESP Simulates localization responses, confidence judgements,
%and unity judgements for a bimodal audiovisual stimulus.
%
% This function models and simulates the outcomes of an experiment where
% subjects are presented with bimodal audiovisual stimuli and are asked to
% make localization responses, judge their confidence in these responses,
% and decide whether the auditory and visual components originated from a
% common source. The simulation is based on parameters defining biases,
% measurement noise, and prior beliefs about stimulus location.
%
% Inputs:
%   pCommon  (double): Common-cause prior probability for vision and audition
%                      being from the same source.
%   num_rep     (int): Number of simulated trials.
%   sA       (double): Location of the auditory component of the AV pair.
%   sV       (double): Location of the visual component of the AV pair.
%   aA       (double): Parameter controlling location-dependent biases in
%                      auditory spatial perception.
%   bA       (double): Parameter controlling location-independent biases in
%                      auditory spatial perception.
%   sigA     (double): Measurement noise for an auditory stimulus.
%   sigV     (double): Measurement noise for a visual stimulus.
%   muP      (double): Center of the prior distribution over stimulus location.
%   sigP     (double): Width of the prior distribution over stimulus location.
%   x        (vector): The screen pixel space, used for simulation scaling or
%                      mapping.
%   fixP     (struct): Fixed parameters for the gain function
%
% Outputs:
%   loc      (matrix): Simulated localization responses. This is a matrix
%                      with dimensions 4 (decision strategies, MA has two 
%                      possibilities: MAP & Expected Value) x 2 (modalities)
%                      x nT (number of trials).
%   conf     (matrix): Simulated confidence judgements for each trial. This
%                      is a matrix with dimensions 4 (decision strategies,
%                       PM has two possibilities) x 2
%                      (modalities) x nT (number of trials).
%   unity    (vector): Simulated unity judgements for each trial. This is a
%                      vector of nT (number of trials).
%
% Date: 24/02 Version: 1.0

x                           = fixP.x;
loc                         = NaN(4,2,num_rep);

%simulate measurements, which are drawn from Gaussian distributions
% stochasticity starts here
mA                          = randn(1, num_rep).*sigA + (sA * aA + bA);
mV                          = randn(1, num_rep).*sigV + sV;

mA                          = bounded(mA,min(x),max(x));
mV                          = bounded(mV,min(x),max(x));

%compute constants (these will come in handy when writing the equations
%for the likelihood/posterior of a common cause and separate causes.
JA                          = sigA^2;
JV                          = sigV^2;
JP                          = sigP^2;
const1                      = JA*JV + JA*JP + JV*JP;
constA                      = JA + JP;
constV                      = JV + JP;

%calculate the likelihood of a common cause and separate causes
%Eq. 4 in Körding et al., 2007
L_C1                        = 1/(2*pi*sqrt(const1)).*exp(-0.5.*((mA - mV).^2 .* ...
    JP + (mV - muP).^2.*JA + (mA - muP).^2.* JV)./const1);
%Eq. 6 in Körding et al., 2007
L_C2                        = 1/(2*pi*sqrt(constA*constV)).*exp(-0.5.*(mA - muP).^2./constA +...
    (mV - muP).^2./constV);

%calculate posterior of a common cause and separate causes
%Eq. 2 in Körding et al., 2007
post_C1                     = pCommon.*L_C1./(pCommon.*L_C1 + (1-pCommon).*L_C2);
%posterior of separate causes = 1 - post_C1
post_C2                     = 1 - post_C1;
% unity judgement
unity                       = post_C1 > 0.5;

%compute the two intermediate location estimates
%An integrated intermediate estimate is the sum of mA, mV and muP with
%each weighted by their relative reliabilities
%Eq. 12 in Körding et al., 2007
sHat_C1                     = (mA./JA + mV./JV + muP/JP)./(1/JV + 1/JA + 1/JP);

% shat_C1 = bounded(shat_C1,min(x),max(x));
%A segregated intermediate estimate is the sum of mA/mV and muP with
%each weighted by their relative reliabilities
%Eq. 11 in Körding et al., 2007
sHat_A_C2                   = (mA./JA + muP/JP)./(1/JA + 1/JP);
sHat_V_C2                   = (mV./JV + muP/JP)./(1/JV + 1/JP);

% initialize intermediate estimate pdfs
pdf.sHat_C1                 = norm_vector(normpdf(repmat(x,num_rep,1),sHat_C1',1/sqrt(1/JV + 1/JA + 1/JP)));
pdf.sHat_A_C2               = norm_vector(normpdf(repmat(x,num_rep,1),sHat_A_C2',1/sqrt(1/JA + 1/JP)));
pdf.sHat_V_C2               = norm_vector(normpdf(repmat(x,num_rep,1),sHat_V_C2',1/sqrt(1/JV + 1/JP)));

%compute the final location estimates if we assume model averaging.
%Based on this strategy, the final location estimate is the sum of the
%two intermediate location estimates, weighted by the corresponding
%causal structure.
%Eq. 4 in Wozny et al., 2010
loc(1,1,:)                  = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
loc(1,2,:)                  = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;

% The posterior for model averaging would be a Gaussian mixture, with its
% expecte value equal to the weighted means of the two causal scenarios.
% But its MAP will not necessarily be there.
pdf.MA_A                    = norm_vector(post_C1'.* pdf.sHat_C1 + post_C2'.* pdf.sHat_A_C2);
pdf.MA_V                    = norm_vector(post_C1'.* pdf.sHat_C1 + post_C2'.* pdf.sHat_V_C2);

[~,loc(2,1,:)]              = max(pdf.MA_A,[],2);
[~,loc(2,2,:)]              = max(pdf.MA_V,[],2);

conf(1,1,:)                 = eGain(pdf.MA_A, round(loc(1,1,:)), fixP.maxPoint, fixP.minPoint, fixP.elbow, fixP.screenX);
conf(1,2,:)                 = eGain(pdf.MA_V, round(loc(1,2,:)), fixP.maxPoint, fixP.minPoint, fixP.elbow, fixP.screenX);
conf(2,1,:)                 = eGain(pdf.MA_A, round(loc(2,1,:)), fixP.maxPoint, fixP.minPoint, fixP.elbow, fixP.screenX);
conf(2,2,:)                 = eGain(pdf.MA_V, round(loc(2,2,:)), fixP.maxPoint, fixP.minPoint, fixP.elbow, fixP.screenX);
%compute the final location estimates if we assume model selection.
%Based on this strategy, the final location estimate depends purely on
%the causal structure that is more probable.
%Eq. 5 in Wozny et al., 2010
loc(3,1:2,:)                = repmat(sHat_C1,[2 1]);
loc(3,1,post_C1<0.5)        = sHat_A_C2(post_C1<0.5);
loc(3,2,post_C1<0.5)        = sHat_V_C2(post_C1<0.5);



% select between intermediate estimated pdfs
slc                         = ~(post_C1<0.5);
pdf.MS_A(slc,:)             = pdf.sHat_C1(slc,:);
pdf.MS_A(~slc,:)            = pdf.sHat_A_C2(~slc,:);
pdf.MS_V(slc,:)             = pdf.sHat_C1(slc,:);
pdf.MS_V(~slc,:)            = pdf.sHat_V_C2(~slc,:);

conf(3,1,:)                 = eGain(pdf.MS_A, round(loc(2,1,:)), fixP.maxPoint, fixP.minPoint, fixP.elbow, fixP.screenX);
conf(3,2,:)                 = eGain(pdf.MS_V, round(loc(2,2,:)), fixP.maxPoint, fixP.minPoint, fixP.elbow, fixP.screenX);

%compute the final location estimates if we assume probability matching.
%Based on this strategy, the final location estimate is the integrated
%one with a probability of post_C1, and is the segregated one with a
%probability of post_C2.
%Eq. 6 in Wozny et al., 2010
idx                         = rand(1,num_rep);
idx_C2                      = (idx > post_C1);
loc(4,1:2,:)                = repmat(sHat_C1,[2 1]);
loc(4,1,idx_C2)             = sHat_A_C2(idx_C2);
loc(4,2,idx_C2)             = sHat_V_C2(idx_C2);

% PM can use the selected intermediate posterior by the posterior of
% common cause
slc2                        = ~idx_C2;
pdf.PM_select_A(slc2,:)     = pdf.sHat_C1(slc2,:);
pdf.PM_select_A(~slc2,:)    = pdf.sHat_A_C2(~slc2,:);
pdf.PM_select_V(slc2,:)     = pdf.sHat_C1(slc2,:);
pdf.PM_select_V(~slc2,:)    = pdf.sHat_V_C2(~slc2,:);
% PM can also use the mixed pdf of intermediate posteriors
pdf.PM_match_A              = norm_vector(pdf.sHat_C1 .* post_C1' + pdf.sHat_A_C2 .* post_C2');
pdf.PM_match_V              = norm_vector(pdf.sHat_C1 .* post_C1' + pdf.sHat_V_C2 .* post_C2');

% conf(3,1,:)                 = eGain(pdf.PM_select_A, round(loc(3,1,:)), fixP.maxPoint, fixP.minPoint, fixP.elbow, fixP.screenX);
% conf(3,2,:)                 = eGain(pdf.PM_select_V, round(loc(3,2,:)), fixP.maxPoint, fixP.minPoint, fixP.elbow, fixP.screenX);
conf(4,1,:)                 = eGain(pdf.PM_match_A, round(loc(3,1,:)), fixP.maxPoint, fixP.minPoint, fixP.elbow, fixP.screenX);
conf(4,2,:)                 = eGain(pdf.PM_match_V, round(loc(3,2,:)), fixP.maxPoint, fixP.minPoint, fixP.elbow, fixP.screenX);

%% unility function
    function v_normed           = norm_vector(v)
        v_normed                    = v ./ sum(v,2);
    end

end