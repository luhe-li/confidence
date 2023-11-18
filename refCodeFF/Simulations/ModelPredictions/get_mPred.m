function [nLL_unity,R] = get_mPred(pC1, epsilon, s_A, s_V, a_A, b_A, sigma_AV_A,...
    sigma_AV_V, CI, mu_P, lapse, model, D_unity)

nLL_unity = 0; 
if isempty(D_unity); bool_fitUnityJdg = 0;
else; bool_fitUnityJdg = 1; end

s_A_prime   = s_A.*a_A + b_A; %the mean of biased auditory measurements
for p = 1:length(s_A_prime)   %for each AV pair with s_A' = s_A_prime(p)
    x1_grid    = linspace(s_A_prime(p) - model.numSD*sigma_AV_A,...
                          s_A_prime(p) + model.numSD*sigma_AV_A,...
                          model.numBins_A);
    mDist_AV_A = norm_dst(x1_grid, s_A_prime(p), sigma_AV_A, 0);
    s_V_prime  = s_V;
    for q = 1:length(s_V_prime) %for each AV pair with s_V' = s_V_prime(q)
        x2_grid    = linspace(s_V_prime(q) - model.numSD*sigma_AV_V,...
                              s_V_prime(q) + model.numSD*sigma_AV_V,...
                              model.numBins_V);
        mDist_AV_V = norm_dst(x2_grid, s_V_prime(q), sigma_AV_V, 0);
                     %normpdf(x2_grid, s_V_prime(q),sigma_AV_V); 
        %use outer product to get the bi-variate distribution
        p_mAmV_given_sAsV = mDist_AV_V'.*mDist_AV_A; 
        p_mAmV_given_sAsV = p_mAmV_given_sAsV./(sum(p_mAmV_given_sAsV(:)));
        %GRID_X1: a matrix where each row is a copy of x1_grid
        %GRID_X2: matrix where each column is a copy of x2_grid
        [GRID_X1, GRID_X2]= meshgrid(x1_grid, x2_grid); 

        %calculate the likelihood and the posterior of a common cause and 
        %separate causes L(C=1|x1=mA, x2=mV) and L(C=2|x1=mA, x2=mV)
        %size of Post_C1 (#rows x #cols): length(x2_grid) x length(x1_grid) 
        [Post_C1, Post_C2, ~, ~]= calculatePostC1C2(GRID_X1, GRID_X2, CI, mu_P, pC1); 
        
        %calculate the shat_{C=1}, shat_{C=2} and MAP estimates assuming
        %model averaging
        shat_C1           = (GRID_X1./CI.J_A + GRID_X2./CI.J_V + ...
                            mu_P./CI.J_P)./CI.constC1_shat;
        [shat_C2, MAP_MA] = deal(NaN(length(model.modality), model.numBins_V,...
                            model.numBins_A));
        shat_C2(1,:,:) = (GRID_X1./CI.J_A + mu_P/CI.J_P)./CI.constC2_1_shat;
        shat_C2(2,:,:) = (GRID_X2./CI.J_V + mu_P/CI.J_P)./CI.constC2_2_shat;
        MAP_MA(1,:,:)  = shat_C1.*Post_C1 + squeeze(shat_C2(1,:,:)).*Post_C2;
        MAP_MA(2,:,:)  = shat_C1.*Post_C1 + squeeze(shat_C2(2,:,:)).*Post_C2;
        
        %---------------------Unity judgment-------------------------------
        if strcmp(model.strategy_unity,'MAP')
            I_C1_resp = squeeze(double(abs(MAP_MA(1,:,:) - MAP_MA(2,:,:)) < epsilon));  
        elseif strcmp(model.strategy_unity, 'measurements')
            I_C1_resp = double(abs(GRID_X1 - GRID_X2) < epsilon);  
        else %'posteriorC1'
            I_C1_resp = double((Post_C1 > 0.5));
        end
        P_C1_resp                 = NaN(size(I_C1_resp));
        P_C1_resp(I_C1_resp == 1) = 1-lapse;
        P_C1_resp(I_C1_resp == 0) = lapse;
        
        if bool_fitUnityJdg
            for k = 1:size(D_unity,3)
                if D_unity(p,q,k) == 1; p_I_given_m = P_C1_resp; 
                else; p_I_given_m = 1 - P_C1_resp; end
                nLL_unity = nLL_unity - log(sum(sum(p_I_given_m.*p_mAmV_given_sAsV)));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%SAVE ONLY IF NEEDED%%%%%%%%%%%%%%%%%%%%%%%
        %save the joint likelihood
        R.p_mAmV_given_sAsV(p,q,:,:) = p_mAmV_given_sAsV;
        %save predictions on the unity judgment
        R.pC1_given_sAsV(p,q) = sum(sum(P_C1_resp.*p_mAmV_given_sAsV));
        %save predictions on localization responses
        if strcmp(model.strategy_MAP,'MA')
            R.MAP(p,q,:,:,:) = MAP_MA;
        elseif strcmp(model.strategy_MAP,'MS') && strcmp(model.strategy_unity,'posteriorC1')
            for m = 1:length(model.modality)
                R.MAP(p,q,m,:,:) = shat_C1;
                R.MAP(p,q,m,Post_C1 < 0.5) = squeeze(shat_C2(m,Post_C1 < 0.5));
            end
        else %measurements
            for m = 1:length(model.modality)
                R.MAP(p,q,m,:,:) = shat_C1;
                R.MAP(p,q,m,abs(GRID_X1 - GRID_X2) > epsilon) = ...
                    squeeze(shat_C2(m,abs(GRID_X1 - GRID_X2) > epsilon));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end


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

function p = norm_dst(x,mu,sigma,t)
p = 1/sqrt(2*pi*sigma^2).*exp(-(x-mu).^2./(2*sigma^2)) + t;
