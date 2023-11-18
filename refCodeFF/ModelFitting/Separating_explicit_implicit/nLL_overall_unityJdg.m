function [nLL, R1, R2] = nLL_overall_unityJdg(a_A, b_A, sigma_A, sigma_V,...
    sigma_AV_A, sigma_AV_V, epsilon, pC1_pre_cong, pC1_post_cong, ...
    pC1_pre_incong, pC1_post_incong, lapse_matching, lapse_unity_cong, ...
    lapse_unity_incong, data, model)
%fixed parameter values
sigma_P = 100; mu_P = 0;
%--------------------------------------------------------------------------
%-------------FOR UNIMODAL LOCALIZATION AND THE MATCHING TASK--------------
%--------------------------------------------------------------------------
%the mean of measurement distributions (auditory measurements are biased)
%data.unimodal: 1st dimension, (1st row: stimulus loc; 2nd row: loc resp)
%               2nd dimension, trial number
%               3rd dimension, modality (1st: A; 2nd: V)
s_A_prime_uni      = a_A.*data.unimodal(1,:,1) + b_A; 
s_V_prime_uni      = data.unimodal(1,:,2);
%to calculate the mean of the response distributions
c_A                = (1/sigma_A^2)/(1/sigma_A^2+1/sigma_P^2);
c_V                = (1/sigma_V^2)/(1/sigma_V^2+1/sigma_P^2);
f_A                = (mu_P/sigma_P^2)/(1/sigma_A^2+1/sigma_P^2);
f_V                = (mu_P/sigma_P^2)/(1/sigma_V^2+1/sigma_P^2);
mu_shat_A_uni      = c_A.*s_A_prime_uni + f_A;
mu_shat_V_uni      = c_V.*s_V_prime_uni + f_V;
R1.mu_shat_A_uni   = unique(mu_shat_A_uni);
R1.mu_shat_V_uni   = unique(mu_shat_V_uni);
%the variances of estimate distributions
sigma_shat_A       = c_A*sigma_A;
sigma_shat_V       = c_V*sigma_V;
R1.sigma_shat_A_wN = sqrt(sigma_shat_A^2 + data.sigma_r^2);
R1.sigma_shat_V_wN = sqrt(sigma_shat_V^2 + data.sigma_r^2);
%unimodal localization task
nLL_unimodal    = calculateNLL_unimodal([mu_shat_A_uni; mu_shat_V_uni], ...
                    [R1.sigma_shat_A_wN; R1.sigma_shat_V_wN], ...
                    [data.unimodal(2,:,1);data.unimodal(2,:,2)]); 
%matching task
%data.matching: 1 by 4 cell, which corresponds to V = -12, -4, 4, 12
%               1st dimension, (1st row: A loc, controlled by staircase)
%                              (2nd row: resp, 1: A right of V; 0: A left)
%               2nd dimension, trial number
nLL_matching    = calculateNLL_matching(a_A, b_A, c_A, c_V, f_A, f_V,...
    data.s_V, sigma_shat_A, sigma_shat_V, lapse_matching, data.matching);

%--------------------------------------------------------------------------
%------------FOR BIMODAL LOCALIZATION (PRE- AND POST-ADAPTATION)-----------
%--------------------------------------------------------------------------
%constants for fitting the unity judgment
CI.J_A            = sigma_AV_A^2;
CI.J_V            = sigma_AV_V^2;
CI.J_P            = sigma_P^2;
CI.constC1        = CI.J_A*CI.J_V + CI.J_A*CI.J_P + CI.J_V*CI.J_P;
CI.constC1_shat   = 1/CI.J_A  + 1/CI.J_V + 1/CI.J_P;
CI.constC2_1      = CI.J_A + CI.J_P;
CI.constC2_1_shat = 1/CI.J_A + 1/CI.J_P;
CI.constC2_2      = CI.J_V + CI.J_P;
CI.constC2_2_shat = 1/CI.J_V + 1/CI.J_P;   

%bimodal localization
%data.bimodal_unity: 1st dim, condition (1: cong; 2: incong)
%                    2nd dim, phase (1: pre; 2: post)
%                    3rd dim, loc A index
%                    4th dim, loc V index
pC1         = [pC1_pre_cong, pC1_post_cong; pC1_pre_incong, pC1_post_incong];
lapse_unity = [lapse_unity_cong, lapse_unity_incong];
R2          = cell(length(model.cond), length(model.phase));
nLL_bimodal = NaN(length(model.cond), length(model.phase));
for i = 1:length(model.cond) 
    lapse_unity_i = lapse_unity(i);
    for j = 1:length(model.phase) 
        Dij_unity = squeeze(data.bimodal_unity(i,j,:,:,:));
        [nLL_bimodal(i,j),R2{i,j}] = calculateNLL_bimodal(pC1(i,j), data.s_A,...
            data.s_V, a_A, b_A, sigma_AV_A, sigma_AV_V, epsilon, CI, mu_P,...
            lapse_unity_i, Dij_unity, data.numUnityTrialsPerLoc, model);
    end
end
%sum the nLL from all tasks together
nLL = nLL_unimodal + nLL_matching + sum(nLL_bimodal(:));


function nLL_unimodal = calculateNLL_unimodal(mu, sig, x)
%we assume that response distributions are Gaussian, centered at mu
%with variance equal to (sigma_shat^2 + sigma_r^2)
%mu, sig, x all consist of 3 rows (1st row: A; 2nd row: V)
LL = arrayfun(@(idx) length(mu(idx,:))*(-0.5*log(2*pi*sig(idx)^2))-...
    sum((x(idx,:) - mu(idx,:)).^2)./(2*sig(idx)^2), 1:2);
%log(1/sqrt(2*pi*sigma^2)*e^(-(x-mu)^2/(2*sigma^2))) =
%-0.5log(2*pi*sigma^2) - (x-mu)^2/(2*sigma^2)
nLL_unimodal = sum(-LL(:));


function nLL_matching = calculateNLL_matching(a_A, b_A, c_A, c_V, f_A,...
    f_V, s_V, sigma_shat_A, sigma_shat_V, lapse, D)
nLL_matching = 0;
%when subjects are presented with a V and an A in order, the probability of
%judging A right of V is 1- \int_{-infty}^0 P(x), where P(x) is a Gaussian,
%centered at (shat_A - shat_V) with variance equal to the sum of both
sig_diff = sqrt(sigma_shat_A^2 + sigma_shat_V^2); 
for j = 1:length(D) %for each standard visual stimulus (there are 4)
    mu_diff = c_A.*(a_A.*D{j}(1,:) + b_A) + f_A - (c_V*s_V(j) + f_V); 
    %probability of reporting A right of V
    pR = lapse/2 + (1-lapse).*(1 - normcdf(0, mu_diff, sig_diff)); 
    LL = sum((D{j}(2,:).*log(pR) + (1-D{j}(2,:)).*log(1-pR)));
    nLL_matching = nLL_matching - LL;
end


function [nLL_bimodal,R] = calculateNLL_bimodal(pC1, s_A, s_V, a_A, b_A, ...
    sigma_AV_A, sigma_AV_V, epsilon, CI, mu_P, lapse, D_unity, N, model)
nLL_bimodal = 0;
s_A_prime   = s_A.*a_A + b_A; %the mean of biased auditory measurements
for p = 1:length(s_A_prime)   %for each AV pair with s_A' = s_A_prime(p)
    x1_grid    = linspace(s_A_prime(p) - model.numSD*sigma_AV_A,...
                          s_A_prime(p) + model.numSD*sigma_AV_A,...
                          model.numBins_A);
    mDist_AV_A = norm_dst(x1_grid, s_A_prime(p), sigma_AV_A, 0);
                 %normpdf(x1_grid, s_A_prime(p), sigma_AV_A);
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
        
        for k = 1:N %N=20
            if D_unity(p,q,k) == 1; p_I_given_m = P_C1_resp; 
            else; p_I_given_m = 1 - P_C1_resp; end
            nLL_bimodal = nLL_bimodal - log(sum(sum(p_I_given_m.*p_mAmV_given_sAsV)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%SAVE ONLY IF NEEDED%%%%%%%%%%%%%%%%%%%%%%%
%         %save the joint likelihood
%         R.p_mAmV_given_sAsV(p,q,:,:) = p_mAmV_given_sAsV;
%         %save predictions on the unity judgment
%         R.pC1_given_sAsV(p,q) = sum(sum(P_C1_resp.*p_mAmV_given_sAsV));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
R = []; %comment it out if we are saving R


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

%figure; imagesc(x1_grid, x2_grid,p_mAmV_given_sAsV);colorbar; 
%xticks(round([x1_grid(1), s_A_prime(p), x1_grid(end)],2)); 
%yticks(round([x2_grid(1),s_V_prime(q), x2_grid(end)],2)); set(gca,'FontSize',15);