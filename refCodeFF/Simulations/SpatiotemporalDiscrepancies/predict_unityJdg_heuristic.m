function  propC1 = predict_unityJdg_heuristic(param, N, spatialD, temporalD,...
    pCommon, ds_unityJdg, c_spatial, c_temporal, c_spatiotemporal, f_s, f_t)

if nargin < 9; c_spatiotemporal = NaN; f_s = NaN; f_t = NaN; end
%useful constants
% CI.Spat_P              = param.sigmaP_spatial^2;
% CI.Spat_V              = param.sigma_spatial_AV_V^2;
% CI.Spat_A              = param.sigma_spatial_AV_A^2;
% CI.Spat_constC1        = CI.Spat_A*CI.Spat_V + CI.Spat_A*CI.Spat_P + ...
%                          CI.Spat_V*CI.Spat_P;
% CI.Spat_constC2_1      = CI.Spat_A + CI.Spat_P;
% CI.Spat_constC2_2      = CI.Spat_V + CI.Spat_P;
% CI.Temp_constC1        = param.sigma_deltaT^2 + param.sigmaP_deltaT_C1^2;
% CI.Temp_constC2        = param.sigma_deltaT^2 + param.sigmaP_deltaT_C2^2;
% CI.Spat_constC1_shat   = 1/CI.Spat_A + 1/CI.Spat_V + 1/CI.Spat_P;
% CI.Spat_constC2_1_shat = 1/CI.Spat_A + 1/CI.Spat_P;
% CI.Spat_constC2_2_shat = 1/CI.Spat_V + 1/CI.Spat_P;

CI.J_A            = param.sigma_spatial_AV_A^2;
CI.J_V            = param.sigma_spatial_AV_V^2;
CI.J_P            = param.sigmaP_spatial^2;
CI.constC1        = CI.J_A*CI.J_V + CI.J_A*CI.J_P + CI.J_V*CI.J_P;
CI.constC1_shat   = 1/CI.J_A  + 1/CI.J_V + 1/CI.J_P;
CI.constC2_1      = CI.J_A + CI.J_P;
CI.constC2_1_shat = 1/CI.J_A + 1/CI.J_P;
CI.constC2_2      = CI.J_V + CI.J_P;
CI.constC2_2_shat = 1/CI.J_V + 1/CI.J_P;   

%initialization
bool_C1  = NaN(1, N);
%for each trial, draw a measurement
m_A_tilde = randn(1, N).* param.sigma_spatial_AV_A;
m_V_tilde = randn(1, N).* param.sigma_spatial_AV_V + spatialD;
m_T       = randn(1, N).* param.sigma_deltaT + temporalD;
for i = 1:N
    %find the posterior of a common cause given those 2 measurements
    [pC1, ~, ~, ~] = calculatePostC1C2(m_A_tilde(i), m_V_tilde(i), CI, param.muP_spatial, pCommon);
    
    %calculate MAP estimates
    [MAP_A, MAP_V] = findMAPestimates(m_A_tilde(i), m_V_tilde(i), pC1, param, CI);
    
    switch ds_unityJdg
        case 'AND'
            if (m_T(i) < c_temporal) && (MAP_A - MAP_V < c_spatial); bool_C1(i) = 1;
            else; bool_C1(i) = 0; end
        case 'OR'
            if (m_T(i) < c_temporal) || (MAP_A - MAP_V < c_spatial); bool_C1(i) = 1;
            else; bool_C1(i) = 0; end
        case 'Euclidean distance'
            if sqrt(f_t*m_T(i)^2 + f_s*(MAP_A - MAP_V)^2) < c_spatiotemporal; bool_C1(i) = 1;
            else; bool_C1(i) = 0; end
    end
end
%no need to return the shifts over time, just the final one is enough.
propC1 = sum(bool_C1)/N;


function [MAP_A, MAP_V] = findMAPestimates(m_A_tilde,m_V_tilde, posteriorC1, param, CI)
MAP_A_C1 = (m_A_tilde/CI.J_A + m_V_tilde/CI.J_V + ...
            param.muP_spatial/CI.J_P)/CI.constC1_shat;
MAP_V_C1 = MAP_A_C1;
MAP_A_C2 = (m_A_tilde/CI.J_A + param.muP_spatial/CI.J_P)/CI.constC2_1_shat;
MAP_V_C2 = (m_V_tilde/CI.J_V + param.muP_spatial/CI.J_P)/CI.constC2_2_shat;
MAP_A    = MAP_A_C1*posteriorC1 + MAP_A_C2*(1-posteriorC1);
MAP_V    = MAP_V_C1*posteriorC1 + MAP_V_C2*(1-posteriorC1);


function [Post_C1, Post_C2, L_C1, L_C2] = calculatePostC1C2(X1, X2, CI, mu_P, pC1)
%likelihood of a common cause and seperate causes
L_C1     = 1/(2*pi*sqrt(CI.constC1))*exp(-0.5*((X1 - X2).^2.*CI.J_P +...
            (X1 - mu_P).^2*CI.J_V + (X2 - mu_P).^2.*CI.J_A)./CI.constC1);
L_C2     = 1/(2*pi*sqrt(CI.constC2_1*CI.constC2_2))*exp(-0.5*...
            ((X1 - mu_P).^2./CI.constC2_1+(X2 - mu_P).^2./CI.constC2_2)); 
normTerm = L_C1.*pC1 + L_C2.*(1-pC1);
%posterior of a common cause
Post_C1  = L_C1.*pC1./normTerm;
Post_C2  = 1 - Post_C1;
