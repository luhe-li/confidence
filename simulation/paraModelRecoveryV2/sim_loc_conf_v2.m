function [out] = sim_loc_conf_v2(sim_model, num_rep, sA, sV,...
    aA, bA, sigA, sigV, sigP, pCommon, sigM, cA, cV,...
    lapse, muP, fixP)

% Date: 24/04/05

x                           = fixP.x;
loc                         = NaN(2,num_rep);

%simulate measurements, which are drawn from Gaussian distributions
% stochasticity starts here
mA                          = randn(1, num_rep).*sigA + (sA * aA + bA);
mV                          = randn(1, num_rep).*sigV + sV;

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
loc(1,:)                  = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
loc(2,:)                  = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;

% calculate posterior variance given each model
switch model_ind
    case 1
        variance(1,:)                 = repmat(1/JA, [1, num_rep]);
        variance(2,:)                 = repmat(1/JV, [1, num_rep]);
    case 2
        variance(1:2,:)               = 1/(1/JV + 1/JA + 1/JP);
        variance(1,post_C1<0.5)       = 1/(1/JA + 1/JP);
        variance(2,post_C1<0.5)       = 1/(1/JV + 1/JP);
    case 3
        variance(1,:)                 = post_C1'.* 1/(1/JV + 1/JA + 1/JP) + post_C2'.* 1/(1/JA + 1/JP) + post_C1'.* post_C2' .* (sHat_A_C2' - sHat_C1').^2;
        variance(2,:)                 = post_C1'.* 1/(1/JV + 1/JA + 1/JP) + post_C2'.* 1/(1/JV + 1/JP) + post_C1'.* post_C2' .* (sHat_V_C2' - sHat_C1').^2;
end


end