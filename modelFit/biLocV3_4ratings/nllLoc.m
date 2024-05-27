function out = nllLoc(freeParam, model, data)

switch model.mode

    case 'initiate'

        out.paraID                   = {'aA','bA','\sigma_{V1}','\sigma_{A}','\sigma_{V2}','\sigma_{P}','p_{common}'};
        out.num_para                 = length(out.paraID);

        % hard bounds, the range for LB, UB, larger than soft bounds
        paraH.aA                     = [ 0.5,     3]; % degree
        paraH.bA                     = [ -10,    10]; % degree
        paraH.sigV1                  = [1e-2,     3]; % degree
        paraH.sigA                   = [ 0.1,    10]; % degree
        paraH.sigV2                  = [ 0.1,    10]; % degree
        paraH.sigP                   = [   1,    30]; % degrees
        paraH.pC1                    = [1e-3,1-1e-3]; % weight

        % soft bounds, the range for PLB, PUB
        paraS.aA                     = [ 0.9,     1]; % degree
        paraS.bA                     = [  -2,     2]; % degree
        paraS.sigV1                  = [ 0.1,     2]; % degree
        paraS.sigA                   = [   3,     8]; % degree
        paraS.sigV2                  = [   3,     8]; % degree
        paraS.sigP                   = [   8,    10]; % degrees
        paraS.pC1                    = [ 0.5,   0.7]; % weight

        % reorganize parameter bounds to feed to bads
        fn                           = fieldnames(paraH);
        for k                        = 1:numel(fn)
            out.lb(:,k)                  = paraH.(fn{k})(1);
            out.ub(:,k)                  = paraH.(fn{k})(2);
            out.plb(:,k)                 = paraS.(fn{k})(1);
            out.pub(:,k)                 = paraS.(fn{k})(2);
        end
        out.paraS                    = paraS; out.paraH = paraH;

        % get grid initializations
        out.init                     = getInit(out.lb, out.ub, model.num_sec, model.num_run);

    case {'optimize','predict'}

        % fixed parameter values for reducing model
        mu_P = 0;
        lapse = 0.02;
        aV = 1;
        bV = 0;

        % free parameters
        aA                           = freeParam(1);
        bA                           = freeParam(2);
        sigV1                        = freeParam(3);
        sigA                         = freeParam(4);
        sigV2                        = freeParam(5);
        sigP                         = freeParam(6);
        pCommon                      = freeParam(7);

        sigVs = [sigV1, sigV2];
        num_sigVs = numel(sigVs);
        sigMotor = data.sigMotor;

        R = cell(1, num_sigVs);
        nLL_bimodal = NaN(1, num_sigVs);

        for ii = 1:num_sigVs

            sigV = sigVs(ii);

            % constants
            CI.J_A            = sigA^2;
            CI.J_V            = sigV^2;
            CI.J_P            = sigP^2;
            CI.constC1        = CI.J_A*CI.J_V + CI.J_A*CI.J_P + CI.J_V*CI.J_P;
            CI.constC1_shat   = 1/CI.J_A  + 1/CI.J_V + 1/CI.J_P;
            CI.constC2_1      = CI.J_A + CI.J_P;
            CI.constC2_1_shat = 1/CI.J_A + 1/CI.J_P;
            CI.constC2_2      = CI.J_V + CI.J_P;
            CI.constC2_2_shat = 1/CI.J_V + 1/CI.J_P;

            % auditory locations (4) x visual locations (4) x postcues (2)
            % x rep
            data_resp = squeeze(data.org_resp(:,:,:,ii,:));

            [nLL_bimodal(ii), R{ii}] = calculateNLLloc(...
                aA, bA, aV, bV, sigA, sigV, pCommon, ...
                CI, mu_P, sigMotor, data_resp, model);

        end

        if strcmp(model.mode, 'optimize')

            out = sum(nLL_bimodal(:));

        elseif strcmp(model.mode, 'predict')

            out = R;

        end

end

end

function [nLL_bimodal, R] = calculateNLLloc(...
    aA, bA, aV, bV, sigA, sigV, pC1,...
    CI, mu_P, sigMotor, data_resp, model)

nLL_bimodal = 0;
sA_prime   = model.sA.*aA + bA; %the mean of biased auditory measurements

for p = 1:length(sA_prime)   %for each AV pair with s_A' = s_A_prime(p)

    x1_grid    = linspace(sA_prime(p) - model.num_SD*sigA,...
        sA_prime(p) + model.num_SD*sigA,...
        model.numBins_A);
    mDist_AV_A = norm_dst(x1_grid, sA_prime(p), sigA, 0);

    sV_prime  = sA_prime; % model.sV.*aV+bV;
    for q = 1:length(sV_prime) %for each AV pair with s_V' = s_V_prime(q)

        x2_grid    = linspace(sV_prime(q) - model.num_SD*sigV,...
            sV_prime(q) + model.num_SD*sigV,...
            model.numBins_V);
        mDist_AV_V = norm_dst(x2_grid, sV_prime(q), sigV, 0);

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

        %----------------------- Compute likelihood -----------------------
        % For each same sA, sV combination, the data are organized by
        % response modality (A, V), repeititon of trials.
        % First half trials are A loc resps, second half trials are V loc
        % resps, from the same sA, sV combination
        locResp_A_V = [squeeze(data_resp(p,q,1,:))'; squeeze(data_resp(p,q,2,:))'];
        [num_mod, num_rep] = size(locResp_A_V);

        if strcmp(model.strategy_loc,'MA') %Model averaging
            for mm = 1:num_mod
                for kk = 1:num_rep
                    p_r_given_MAP = norm_dst(locResp_A_V(mm, kk),squeeze(MAP_MA(mm,:,:)),...
                        sigMotor,1e-20);
                    nLL_bimodal = nLL_bimodal - log(sum(sum(p_r_given_MAP.*p_mAmV_given_sAsV)));
                end
            end
        elseif strcmp(model.strategy_loc,'MS') %Model selection
            for mm = 1:num_mod
                for kk = 1:num_rep
                    if Post_C1(mm,kk) > 0.5
                        p_r_given_MAP = norm_dst(locResp_A_V(mm, kk),shat_C1,sigMotor,1e-20);
                    else
                        p_r_given_MAP = norm_dst(locResp_A_V(mm, kk),squeeze(shat_C2(mm,:,:)),sigMotor,1e-20);
                    end
                    nLL_bimodal = nLL_bimodal - log(sum(sum(p_r_given_MAP.*p_mAmV_given_sAsV)));
                end
            end

        end

        %----------------------- Save if predicted -----------------------
        R = [];
        if strcmp(model.mode,'predict')
            
            %save the joint likelihood
            R.p_mAmV_given_sAsV(p,q,:,:) = p_mAmV_given_sAsV;
            % save the probability of confidence report
            for mm = 1:length(model.modality)
                R.p_conf_given_sAsV(p,q,mm) = sum(sum(squeeze(temp_p_conf(mm,:,:)).*p_mAmV_given_sAsV));
            end

            %save predictions on location estimates
            if strcmp(model.strategy_loc,'MA')
                R.loc(p,q,:,:,:) = MAP_MA;
            elseif strcmp(model.strategy_loc,'MS')
                for mm = 1:length(model.modality)
                    R.loc(p,q,mm,:,:) = shat_C1;
                    R.loc(p,q,mm,Post_C1 < 0.5) = squeeze(shat_C2(mm,Post_C1 < 0.5));
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
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
end

function p = norm_dst(x,mu,sigma,t)
p = 1/sqrt(2*pi*sigma^2).*exp(-(x-mu).^2./(2*sigma^2)) + t;
end