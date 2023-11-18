function  [propC1, pC1_end, pC1_update] = simulate_learningPhase(...
    param, alpha_pC1, D, model)
%initialize shifts in the past and set them to 0 when t = 1 
pC1_update     = zeros(1,D.totalTrials+1);
pC1_update(1)  = param.pCommon;
bool_C1        = NaN(1, D.totalTrials);

%constants for fitting the unity judgment
CI.J_A            = param.sigma_A_learning^2;
CI.J_V            = param.sigma_V_learning^2;
CI.J_P            = param.sigma_P^2;
CI.constC1        = CI.J_A*CI.J_V + CI.J_A*CI.J_P + CI.J_V*CI.J_P;
CI.constC1_shat   = 1/CI.J_A  + 1/CI.J_V + 1/CI.J_P;
CI.constC2_1      = CI.J_A + CI.J_P;
CI.constC2_1_shat = 1/CI.J_A + 1/CI.J_P;
CI.constC2_2      = CI.J_V + CI.J_P;
CI.constC2_2_shat = 1/CI.J_V + 1/CI.J_P; 

for i = 1:D.totalTrials    
    %for each trial, draw a measurement
    m_A_tilde = randn * param.sigma_A_learning + (param.a_A*D.AVpairs(1,i) + param.b_A);
    m_V_tilde = randn * param.sigma_V_learning + D.AVpairs(2,i);

    %find the posterior of a common cause given those 2 measurements
    [Post_C1, ~, ~, ~] = calculatePostC1C2(m_A_tilde, m_V_tilde, CI,...
        pC1_update(i), param);
    
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
                if Post_C1 > .5; bool_C1(i) = 1; 
                else; bool_C1(i) = 0; end 
            case 'MAP'
                if abs(MAP_A - MAP_V) < param.c_unityJdg; bool_C1(i) = 1;
                else; bool_C1(i) = 0; end
        end
    end
    pC1_update(i+1)  = pC1_update(i) + alpha_pC1*(Post_C1 - pC1_update(i));
end

%no need to return the shifts over time, just the final one is enough.
propC1  = nansum(bool_C1)/sum(~isnan(bool_C1));
pC1_end = pC1_update(end);


function [Post_C1, Post_C2, L_C1, L_C2] = calculatePostC1C2(X1, X2, CI, pC1, param)
%likelihood of a common cause and seperate causes
L_C1     = 1/(2*pi*sqrt(CI.constC1))*exp(-0.5*((X1 - X2).^2.*CI.J_P +...
            (X1 - param.mu_P).^2*CI.J_V + (X2 - param.mu_P).^2.*CI.J_A)./CI.constC1);
L_C2     = 1/(2*pi*sqrt(CI.constC2_1*CI.constC2_2))*exp(-0.5*...
            ((X1 - param.mu_P).^2./CI.constC2_1+(X2 - param.mu_P).^2./CI.constC2_2)); 
normTerm = L_C1.*pC1 + L_C2.*(1-pC1);
%posterior of a common cause
Post_C1  = L_C1.*pC1./normTerm;
Post_C2  = 1 - Post_C1;


function [MAP_A, MAP_V] = findMAPestimates(m_A_tilde, m_V_tilde, pC1,...
    postC1, CI, param, model)

MAP_C1   = (m_A_tilde/CI.J_A + m_V_tilde/CI.J_V +  param.mu_P/CI.J_P)./CI.constC1_shat;
MAP_A_C2 = (m_A_tilde/CI.J_A + param.mu_P/CI.J_P)/CI.constC2_1_shat;
MAP_V_C2 = (m_V_tilde./CI.J_V + param.mu_P/CI.J_P)./CI.constC2_2_shat;
                        
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
