function  [propC1, pC1_end,pC1_update] = simulate_learningPhase_wTiming(...
    param, alpha_pC1, D, model)
%initialize shifts in the past and set them to 0 when t = 1
pC1_update      = zeros(1,D.totalTrials+1); %1: spatial; 2: temporal
pC1_update(1)   = param.pCommon;
bool_C1         = NaN(1, D.totalTrials);

%useful constants
CI.Spat_A              = param.sigma_A_learning^2;
CI.Spat_V              = param.sigma_V_learning^2;
CI.Spat_P              = param.sigma_P^2;
CI.Spat_constC1        = CI.Spat_A*CI.Spat_V + CI.Spat_A*CI.Spat_P + ...
                         CI.Spat_V*CI.Spat_P;
CI.Spat_constC2_1      = CI.Spat_A + CI.Spat_P;
CI.Spat_constC2_2      = CI.Spat_V + CI.Spat_P;
CI.Temp_constC1        = param.sigma_deltaT^2 + param.sigmaP_deltaT_C1^2;
CI.Temp_constC2        = param.sigma_deltaT^2 + param.sigmaP_deltaT_C2^2;
CI.Spat_constC1_shat   = 1/CI.Spat_A + 1/CI.Spat_V + 1/CI.Spat_P;
CI.Spat_constC2_1_shat = 1/CI.Spat_A + 1/CI.Spat_P;
CI.Spat_constC2_2_shat = 1/CI.Spat_V + 1/CI.Spat_P;

for i = 1:D.totalTrials
    %for each trial, draw a measurement
    m_A_tilde = randn * param.sigma_A_learning + (param.a_A*D.AVpairs(1,i) + param.b_A);
    m_V_tilde = randn * param.sigma_V_learning + D.AVpairs(2,i);
    m_T       = randn * param.sigma_deltaT + D.AVpairs(4,i) - D.AVpairs(3,i);

    %find the posterior of a common cause given those 2 measurements
    [Post_C1, ~, ~, ~] = findPosteriorC1_spatiotemporal(pC1_update(i),m_A_tilde,...
        m_V_tilde, m_T, param, CI);
    
    %calculate MAP estimates
    [MAP_A, MAP_V] = findMAPestimates(m_A_tilde, m_V_tilde, pC1_update(i),...
        Post_C1, CI, param, model);

    %update the shifts
    if D.bool_unity(i)
        switch model.M_unityJdg
            case 'measurements'
                if abs(m_A_tilde - m_V_tilde) < param.c_unityJdg; bool_C1(i) = 1; 
                else; bool_C1(i) = 0; end
            case 'posteriorC1'
                if Post_C1 > param.c_unityJdg; bool_C1(i) = 1; 
                else; bool_C1(i) = 0; end 
            case 'MAP'
                if abs(MAP_A - MAP_V) < param.c_unityJdg; bool_C1(i) = 1;
                else; bool_C1(i) = 0; end
        end
    end
    pC1_update(i+1) = pC1_update(i) + alpha_pC1*(Post_C1 - pC1_update(i));
end
%no need to return the shifts over time, just the final one is enough.
propC1  = nansum(bool_C1)/sum(~isnan(bool_C1));
pC1_end = pC1_update(end);

function [Post_C1, Post_C2,Spat_L_C1,Spat_L_C2] = findPosteriorC1_spatiotemporal(...
    pC1, m_A_tilde, m_V_tilde, m_T, param,CI)
%Calculate the likelihood of C=1 given m_A and m_V  
Spat_L_C1   = 1/(2*pi*sqrt(CI.Spat_constC1))*exp(-0.5*(...
                (m_A_tilde - m_V_tilde)^2*CI.Spat_P +...
                (m_A_tilde - param.mu_P)^2*CI.Spat_V + ...
                (m_V_tilde - param.mu_P)^2*CI.Spat_A)/CI.Spat_constC1);

Spat_L_C2   = 1/(2*pi*sqrt(CI.Spat_constC2_1*CI.Spat_constC2_2))*...
                exp(-0.5*((m_A_tilde - param.mu_P)^2/CI.Spat_constC2_1+...
                (m_V_tilde - param.mu_P)^2/CI.Spat_constC2_2));

Temp_L_C1   = 1/sqrt(2*pi*CI.Temp_constC1)*exp(-0.5*(m_T - ...
                param.mu_P_deltaT_C1)^2/CI.Temp_constC1);
Temp_L_C2   = 1/sqrt(2*pi*CI.Temp_constC2)*exp(-0.5*(m_T - ...
                param.mu_P_deltaT_C2)^2/CI.Temp_constC2);

%calculate the posterior of C=1 and C=2
normTerm    = Spat_L_C1*Temp_L_C1*pC1 + Spat_L_C2*Temp_L_C2*(1-pC1);
Post_C1     = Spat_L_C1*Temp_L_C1*pC1./normTerm;
Post_C2     = 1 - Post_C1;

function [MAP_A, MAP_V] = findMAPestimates(m_A_tilde, m_V_tilde, pC1,...
    postC1, CI, param, model)

MAP_C1   = (m_A_tilde/CI.Spat_A + m_V_tilde/CI.Spat_V +  param.mu_P/CI.Spat_P)./CI.Spat_constC1_shat;
MAP_A_C2 = (m_A_tilde/CI.Spat_A + param.mu_P/CI.Spat_P)/CI.Spat_constC2_1_shat;
MAP_V_C2 = (m_V_tilde./CI.Spat_V + param.mu_P/CI.Spat_P)./CI.Spat_constC2_2_shat;
                        
switch model.M_estimates
    case 'MA'
        MAP_A = MAP_C1*pC1 + MAP_A_C2*(1-pC1);
        MAP_V = MAP_C1*pC1 + MAP_V_C2*(1-pC1);
    case 'MS'
        if strcmp(model.M_unityJdg, 'measurements')==1
            if abs(m_A_tilde - m_V_tilde) < model.M_unityJdg
                MAP_A = MAP_C1; MAP_V = MAP_C1;
            else; MAP_A = MAP_A_C2; MAP_V = MAP_V_C2;
            end
        elseif strcmp(model.M_unityJdg, 'posteriorC1')==1
            if postC1 > 0.5; MAP_A = MAP_C1; MAP_V = MAP_C1;
            else; MAP_A = MAP_A_C2; MAP_V = MAP_V_C2;end
        end
end
