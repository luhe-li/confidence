function [loc, conf, variance] = sim_loc_pconf(num_rep, sA, sV, aA, bA, sigA, sigV, sigP, pCommon, criterion, fixP)

%sim_loc_pconf simulates localization responses, confidence judgements,
%for a bimodal audiovisual stimulus.
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
%   criterion(double): Confidence criterion on variance of the posterior
%                      distribution
%
% Outputs:
%   loc      (matrix): Simulated localization responses. This is a matrix
%                      with dimensions 4 (decision strategies, MA has two 
%                      possibilities: MAP & Expected Value) x 2 (modalities)
%                      x nT (number of trials).
%   conf     (matrix): Simulated confidence judgements for each trial. 
%                      All using EV loc as estimate location. This
%                      is a matrix with dimensions 3 (three models, 
%                      m1 = using unisensory uncertainty, 
%                      m2 = using a selected intermediate uncertainty, 
%                      m3 = using the full mixture posterior) 
%                      x 2 (modalities) x nT (number of trials).
%
% Date: 24/03 Version: 3.0
% norm_vector = @(v) v./sum(v,2);
f_logit = @(x) 1 ./ (1 + exp(-x));
bounded = @(x, LB, UB) max(LB,min(x, UB)); 

x                           = fixP.x;
loc                         = NaN(3,2,num_rep);
muP = 0;

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

%compute the final location estimates if we assume model averaging.
%Based on this strategy, the final location estimate is the sum of the
%two intermediate location estimates, weighted by the corresponding
%causal structure.
%Eq. 4 in Wozny et al., 2010
loc(1,1,:)                  = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
loc(1,2,:)                  = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;

% calculate posterior variance given each model
variance(1,1,:)                 = post_C1'.* 1/(1/JA + 1/JP) + post_C2'.* 1/(1/JV + 1/JA + 1/JP) + post_C1'.* post_C2' .* (sHat_A_C2' - sHat_C1').^2;
variance(1,2,:)                 = post_C1'.* 1/(1/JV + 1/JP) + post_C2'.* 1/(1/JV + 1/JA + 1/JP) + post_C1'.* post_C2' .* (sHat_V_C2' - sHat_C1').^2;

variance(2,1:2,:)                = 1/(1/JV + 1/JA + 1/JP);
variance(2,1,post_C1<0.5)        = 1/(1/JA + 1/JP);
variance(2,2,post_C1<0.5)        = 1/(1/JV + 1/JP);

variance(3,1,:)                 = JA;
variance(3,2,:)                 = JV;

% conf_var = f_logit(variance);
conf = variance > repmat(criterion, [1, 2, num_rep]);%  model, 2 modalities, num_rep

end