function  [propC1, pCommon_end,pCommon_update] = simulateLearningPhase_update_pC1_maskingEffects(...
    param, D, cond, ds_unityJdg, c_unityJdg, ds_locResp)
%initialize shifts in the past and set them to 0 when t = 1
pCommon_update      = zeros(1,D.totalTrials+1); %1: spatial; 2: temporal
pCommon_update(1)   = param.pCommon;
sigma_spatial_AV_V  = param.sigma_spatial_AV_V;
sigma_spatial_AV_A  = param.sigma_spatial_AV_A;
bool_C1             = NaN(1, D.totalTrials);
for i = 1:D.totalTrials
    %for each trial, draw a measurement
    A_time_relative = D.AVpairs(cond,2,i);
    V_time_relative = D.AVpairs(cond,4,i);
    m_T             = randn * param.sigma_deltaT + A_time_relative - V_time_relative;
    if cond == 2 %incongruent condition
        if m_T < 0 %A precedes V
            param.sigma_spatial_AV_V = sigma_spatial_AV_V*(1+param.gamma_AV_V);
            param.sigma_spatial_AV_A = sigma_spatial_AV_A;
        else
            param.sigma_spatial_AV_A = sigma_spatial_AV_A*(1+param.gamma_AV_A);
            param.sigma_spatial_AV_V = sigma_spatial_AV_V;
        end
    else %congruent condition
        param.sigma_spatial_AV_V = max([sigma_spatial_AV_V*(1-param.gamma_AV_V),0.01]);
        param.sigma_spatial_AV_A = max([sigma_spatial_AV_A*(1-param.gamma_AV_A),0.01]);
    end
    m_A_tilde = randn * param.sigma_spatial_AV_A + (param.a_A*D.AVpairs(cond,1,i) + param.b_A);
    m_V_tilde = randn * param.sigma_spatial_AV_V + D.AVpairs(cond,3,i);

    %find the posterior of a common cause given those 2 measurements
    [pC1, ~, ~, ~] = findPosteriorC1_spatiotemporal(...
        pCommon_update(i),m_A_tilde, m_V_tilde, m_T, param);
    
    %calculate MAP estimates
    [MAP_A, MAP_V] = findMAPestimates(ds_locResp, ds_unityJdg,c_unityJdg, ...
        m_A_tilde, m_V_tilde,pCommon_update(i),pC1, param);

    %update the shifts
    switch ds_unityJdg
        case 'measurements'
            if abs(m_A_tilde - m_V_tilde) < c_unityJdg; bool_C1(i) = 1; 
            else; bool_C1(i) = 0; end
        case 'posteriorC1'
            if pC1 > c_unityJdg; bool_C1(i) = 1; 
            else; bool_C1(i) = 0; end 
        case 'MAP'
            if abs(MAP_A - MAP_V) < c_unityJdg; bool_C1(i) = 1;
            else; bool_C1(i) = 0; end
    end

    pCommon_update(i+1) = pCommon_update(i) + param.alpha*(pC1 - pCommon_update(i));
end
%no need to return the shifts over time, just the final one is enough.
propC1 = sum(bool_C1)/D.totalTrials;
pCommon_end = pCommon_update(end);

function [posteriorC1, posteriorC2,Spat_L_C1,Spat_L_C2] = findPosteriorC1_spatiotemporal(...
    pCommon, m_A_tilde, m_V_tilde, m_T, param)
%Calculate the likelihood of C=1 given m_A and m_V  
Spat_P         = param.sigmaP_spatial^2;
Spat_V         = param.sigma_spatial_AV_V^2;
Spat_A         = param.sigma_spatial_AV_A^2;
Spat_constC1   = Spat_A*Spat_V + Spat_A*Spat_P + Spat_V*Spat_P;
Spat_constC2_1 = Spat_A + Spat_P;
Spat_constC2_2 = Spat_V + Spat_P;
Temp_constC1   = param.sigma_deltaT^2 + param.sigmaP_deltaT_C1^2;
Temp_constC2   = param.sigma_deltaT^2 + param.sigmaP_deltaT_C2^2;

Spat_L_C1   = 1/(2*pi*sqrt(Spat_constC1))*exp(-0.5*(...
                (m_A_tilde - m_V_tilde)^2*Spat_P +...
                (m_A_tilde - param.muP_spatial)^2*Spat_V + ...
                (m_V_tilde - param.muP_spatial)^2*Spat_A)/Spat_constC1);

Spat_L_C2   = 1/(2*pi*sqrt(Spat_constC2_1*Spat_constC2_2))*...
                exp(-0.5*((m_A_tilde - param.muP_spatial)^2/Spat_constC2_1+...
                (m_V_tilde - param.muP_spatial)^2/Spat_constC2_2));

Temp_L_C1   = 1/sqrt(2*pi*Temp_constC1)*exp(-0.5*(m_T - param.muP_deltaT_C1)^2/Temp_constC1);
Temp_L_C2   = 1/sqrt(2*pi*Temp_constC2)*exp(-0.5*(m_T - param.muP_deltaT_C2)^2/Temp_constC2);

%calculate the posterior of C=1 and C=2
posteriorC1 = Spat_L_C1*Temp_L_C1*pCommon/(Spat_L_C1*Temp_L_C1*pCommon + ...
                Spat_L_C2*Temp_L_C2*(1-pCommon));
posteriorC2 = 1 - posteriorC1;

function [MAP_A, MAP_V] = findMAPestimates(ds_locResp, ds_unityJdg,c_unityJdg, m_A_tilde,...
    m_V_tilde,pCommon,posteriorC1, param)
Spat_P         = param.sigmaP_spatial^2;
Spat_V         = param.sigma_spatial_AV_V^2;
Spat_A         = param.sigma_spatial_AV_A^2;
Spat_constC1_shat   = 1/Spat_A + 1/Spat_V + 1/Spat_P;
Spat_constC2_1_shat = 1/Spat_A + 1/Spat_P;
Spat_constC2_2_shat = 1/Spat_V + 1/Spat_P;   

MAP_A_C1 = (m_A_tilde/Spat_A + m_V_tilde/Spat_V + ...
                    param.muP_spatial/Spat_P)/Spat_constC1_shat;
MAP_V_C1 = MAP_A_C1;
MAP_A_C2 = (m_A_tilde/Spat_A + param.muP_spatial/Spat_P)/Spat_constC2_1_shat;
MAP_V_C2 = (m_V_tilde/Spat_V + param.muP_spatial/Spat_P)/Spat_constC2_2_shat;
switch ds_locResp
    case 'MA'
        MAP_A = MAP_A_C1*pCommon + MAP_A_C2*(1-pCommon);
        MAP_V = MAP_V_C1*pCommon + MAP_V_C2*(1-pCommon);
    case 'MS'
        if strcmp(ds_unityJdg, 'measurements')==1
            if abs(m_A_tilde - m_V_tilde) < c_unityJdg
                MAP_A = MAP_A_C1; MAP_V = MAP_V_C1;
            else; MAP_A = MAP_A_C2; MAP_V = MAP_V_C2;end
        elseif strcmp(ds_unityJdg, 'posteriorC1')==1
            if posteriorC1 > 0.5; MAP_A = MAP_A_C1; MAP_V = MAP_V_C1;
            else; MAP_A = MAP_A_C2; MAP_V = MAP_V_C2;end
        end
end

%figure; plot(log(spatialLC1)); hold on; plot(log(spatialLC2)); 
%legend({'log(P(m_A,m_V|C=1))','log(P(m_A,m_V|C=2))'}); box off; set(gca,'FontSize',15)
