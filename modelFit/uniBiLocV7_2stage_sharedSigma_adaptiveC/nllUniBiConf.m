function out = nllUniBiConf(freeParam, model, data)

switch model.mode

    case 'initiate'

        out.paraID                   = {'\sigma_{C}','c1','\delta_{c2}','\delta_{c3}'};
        out.num_para                 = length(out.paraID);

        % hard bounds, the range for LB, UB, larger than soft bounds
        paraH.sigC                   = [1e-4,     5]; % measurement noise of confidence
        paraH.c1                     = [1e-3,     5];
        paraH.dc2                    = [0.01,     5];
        paraH.dc3                    = [0.01,     5];

        % soft bounds, the range for PLB, PUB
        paraS.sigC                   = [1e-4,     2]; % measurement noise of confidence
        paraS.c1                     = [   1,     2];
        paraS.dc2                    = [ 0.1,   0.5];
        paraS.dc3                    = [ 0.1,   0.5];
        paraS.c4                     = [   1,     2];
        paraS.dc5                    = [ 0.1,   0.5];
        paraS.dc6                    = [ 0.1,   0.5];

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

        % free parameters
        sigC                         = freeParam(1);
        c1                           = freeParam(2);
        dc2                          = freeParam(3);
        dc3                          = freeParam(4);
        c4                           = freeParam(5);
        dc5                          = freeParam(6);
        dc6                          = freeParam(7);

        % fixed parameters from first-stage fit
        aA                           = model.locP(1);
        bA                           = model.locP(2);
        sigV1                        = model.locP(3);
        sigA                         = model.locP(4);
        sigV2                        = model.locP(5);
        sigP                         = model.locP(6);
        pCommon                      = model.locP(7);

        % convert
        sigVs = [sigV1, sigV2];
        num_sigVs = numel(sigVs);
        c2 = c1 + dc2;
        c3 = c1 + dc2 + dc3;
        c5 = c4 + dc5;
        c6 = c4 + dc5 + dc6;

        % freeze parameters from data
        sigMotor = data.sigMotor;
        lapse = 0.06;
        aV = 1;
        bV = 0;
        muP = 0;

        if strcmp(model.mode, 'optimize')

            %% unimodal confidence

            % normalize variance by the corresponding modality noise to
            % approximate uncertainty
            norm_var_A = 1/(1/sigP^2 + 1/sigA^2)/sigA;
            norm_var_V1 = 1/(1/sigP^2 + 1/sigV1^2)/sigV1;
            norm_var_V2 = 1/(1/sigP^2 + 1/sigV2^2)/sigV2;

            % unimodal confidence task
            nLL_uni_conf = calculateNLL_uniConf([norm_var_A, norm_var_V1, norm_var_V2],...
                sigC, [c1, c4, c4], [c2, c5, c5], [c3, c6, c6], lapse, data.uni_conf);

            %% bimodal confidence

            nLL_bimodal = NaN(1, num_sigVs);

            for vv = 1:num_sigVs

                sigV = sigVs(vv);

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
                data_resp = squeeze(data.bi_loc(:,:,:,vv,:));
                data_conf = squeeze(data.bi_conf(:,:,:,vv,:));

                [nLL_bimodal(vv)] = calculateNLL_bimodal(...
                    aA, bA, aV, bV, sigA, sigV, sigC, pCommon, [c1, c4], [c2, c5], [c3, c6],...
                    CI, muP, sigMotor, lapse, data_resp, data_conf, model);

            end

            %% sum nll

            out = nLL_uni_conf + sum(nLL_bimodal(:));

        elseif strcmp(model.mode, 'predict')

            fixP.uni_sA = model.uni_sA;
            fixP.uni_sV = model.uni_sV;
            comb_s = combvec(model.bi_sA, model.bi_sV);
            fixP.bi_sA = comb_s(1,:);
            fixP.bi_sV = comb_s(2,:);
            fixP.uni_nrep = model.uni_nrep;
            fixP.bi_nrep = model.bi_nrep;
            fixP.model_ind = model.model_slc;
            fixP.sigMotor = sigMotor;

            [out.uni_loc, out.uni_conf, out.uni_est_var] = simUni(...
                aA, bA, sigA, sigV1, sigV2, muP, sigP, sigC, [c1; c4; c4], [c2; c5; c5], [c3; c6; c6], lapse, fixP);

            [out.bi_loc, out.bi_conf, ~, ~, out.bi_est_var, out.criteria] = simAllModels5D(...
                aA, bA, sigA, [sigV1, sigV2], muP, sigP, pCommon, sigC, c1, c2, c3, lapse, fixP);

            %             [bi_loc, bi_conf] = deal(NaN(numel(bi_sA), numel(bi_sV), numel(model.modality), numel(sigVs), model.bi_nrep));
% 
% 
% 
%             for aa = 1:numel(bi_sA)
%                 for vv = 1:numel(bi_sV)
%                     for rr = 1:numel(sigVs)
% 
%                         fixP.bi_sA = bi_sA(aa);
%                         fixP.bi_sV = bi_sV(vv);
% 
%                         [bi_loc(aa,vv,:,rr,:), bi_conf(aa,vv,:,rr,:), ~, ~, bi_est_var(aa, vv, :, rr, :)] = simAllModels(...
%                             aA, bA, sigA, sigVs(rr), muP, sigP, pCommon, sigC, [c1; c4], [c2; c5], [c3; c6], lapse, fixP);
% 
%                     end
%                 end
%             end
% 
%             out.bi_loc = bi_loc;
%             out.bi_conf = bi_conf;
%             out.bi_est_var = bi_est_var;

        end
end


    function nLL_uni_conf = calculateNLL_uniConf(m, v, c1, c2, c3, lapse, data)

        % convert to lognormal parameters
        mu = log((m.^2)./sqrt(v+m.^2))';
        sigma = sqrt(log(v./(m.^2)+1))';
        temp_p4 = logncdf(c1', mu, sigma);
        temp_p3 = logncdf(c2', mu, sigma) - logncdf(c1', mu, sigma);
        temp_p2 = logncdf(c3', mu, sigma) - logncdf(c2', mu, sigma);

        % add lapse
        p4_conf = lapse./4 + (1-lapse) .* temp_p4 + 1e-20;
        p3_conf = lapse./4 + (1-lapse) .* temp_p3 + 1e-20;
        p2_conf = lapse./4 + (1-lapse) .* temp_p2 + 1e-20;
        p1_conf = 1 - p2_conf - p3_conf - p4_conf;

        % organize probability into a matrix index (A, V1, V2 x p1, p2, p3, p4)
        p_conf = [p1_conf, p2_conf, p3_conf, p4_conf];

        for modality = 1:3
            p = p_conf(modality,:);
            LL(modality, :) = log(p(data(modality,:)));
        end

        nLL_uni_conf = sum(-LL(:));
    end


    function [nLL_bimodal] = calculateNLL_bimodal(...
            aA, bA, aV, bV, sigA, sigV, sigC, pCommon, c1, c2, c3,...
            CI, mu_P, sigMotor, lapse, data_resp, data_conf, model)

        nLL_bimodal = 0;
        sA_prime   = model.bi_sA.*aA + bA; %the mean of biased auditory measurements

        for p = 1:length(sA_prime)   %for each AV pair with s_A' = s_A_prime(p)

            x1_grid    = linspace(sA_prime(p) - model.num_SD*sigA,...
                sA_prime(p) + model.num_SD*sigA,...
                model.numBins_A);
            mDist_AV_A = norm_dst(x1_grid, sA_prime(p), sigA, 0);

            sV_prime  = model.bi_sV.*aV+bV;
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
                [Post_C1, Post_C2, ~, ~]= calculatePostC1C2(GRID_X1, GRID_X2, CI, mu_P, pCommon);

                %calculate the shat_{C=1}, shat_{C=2} and MAP estimates assuming
                %model averaging
                shat_C1           = (GRID_X1./CI.J_A + GRID_X2./CI.J_V + ...
                    mu_P./CI.J_P)./CI.constC1_shat;
                [shat_C2, ~] = deal(NaN(length(model.modality), model.numBins_V,...
                    model.numBins_A));
                shat_C2(1,:,:) = (GRID_X1./CI.J_A + mu_P/CI.J_P)./CI.constC2_1_shat;
                shat_C2(2,:,:) = (GRID_X2./CI.J_V + mu_P/CI.J_P)./CI.constC2_2_shat;
                %                 MAP_MA(1,:,:)  = shat_C1.*Post_C1 + squeeze(shat_C2(1,:,:)).*Post_C2;
                %                 MAP_MA(2,:,:)  = shat_C1.*Post_C1 + squeeze(shat_C2(2,:,:)).*Post_C2;

                %--------------------- Confidence report --------------------------
                % three strategies differ in the variance that is used to derive
                % the confidence variable
                if strcmp(model.strategy_conf,'Optimal')

                    var(1,:,:) = Post_C1 ./CI.constC1_shat + ...
                        Post_C2 ./CI.constC2_1_shat + ...
                        Post_C1 .* Post_C2 .* (squeeze(shat_C2(1,:,:)) - shat_C1).^2;
                    var(2,:,:) = Post_C1 ./CI.constC1_shat + ...
                        Post_C2 ./CI.constC2_2_shat + ...
                        Post_C1 .* Post_C2 .* (squeeze(shat_C2(2,:,:)) - shat_C1).^2;

                elseif strcmp(model.strategy_conf, 'Suboptimal')

                    var(1,:,:)= Post_C1 ./CI.constC1_shat + Post_C2 ./CI.constC2_1_shat;
                    var(2,:,:)= Post_C1 ./CI.constC1_shat + Post_C2 ./CI.constC2_2_shat;
                    %             var = repmat(1/CI.constC1_shat, [2, size(Post_C1)]);
                    %             indices = (Post_C1 >= 0.5);
                    %             var(1,indices) = 1/CI.constC2_1_shat;
                    %             var(2,indices) = 1/CI.constC2_2_shat;

                elseif strcmp(model.strategy_conf, 'Heuristic')

                    var(1,:,:) = repmat(CI.J_P * sigA/CI.constC2_1, size(Post_C1));
                    var(2,:,:) = repmat(CI.J_P * sigV/CI.constC2_2, size(Post_C1));
                    %  var(1,:,:) = repmat(CI.J_A, size(Post_C1));
                    %  var(2,:,:) = repmat(CI.J_V, size(Post_C1));

                end

                % normalize variance by the corresponding modality noise
                norm_var(1,:,:) = var(1,:,:)./sigA;
                norm_var(2,:,:) = var(2,:,:)./sigV;

                % probability of reporting confidence is the value of a lognormal
                % cumulative distribution with a mean of confidence variable (var)
                % and an S.D. of confidence measurement noise (sigC), evaluated at
                % each criteria.
                m = norm_var;
                v = sigC;
                mu = log((m.^2)./sqrt(v+m.^2));
                sigma = sqrt(log(v./(m.^2)+1));
                temp_p4 = logncdf(c1', mu, sigma);
                temp_p3 = logncdf(c2', mu, sigma) - logncdf(c1', mu, sigma);
                temp_p2 = logncdf(c3', mu, sigma) - logncdf(c2', mu, sigma);

                % add lapse
                p4_conf = lapse./4 + (1-lapse) .* temp_p4 + 1e-20;
                p3_conf = lapse./4 + (1-lapse) .* temp_p3 + 1e-20;
                p2_conf = lapse./4 + (1-lapse) .* temp_p2 + 1e-20;
                p1_conf = 1 - p2_conf - p3_conf - p4_conf;

                %----------------------- Compute likelihood -----------------------
                % For each same sA, sV combination, the data are organized by
                % response modality (A, V), repeititon of trials.
                % First half trials are A loc resps, second half trials are V loc
                % resps, from the same sA, sV combination
                confResp_A_V = [squeeze(data_conf(p,q,1,:))'; squeeze(data_conf(p,q,2,:))'];
                [num_modality, num_rep] = size(confResp_A_V);

                if strcmp(model.strategy_loc,'MA') %Model averaging
                    for mm = 1:num_modality
                        for kk = 1:num_rep

                            % confidence probability
                            if confResp_A_V(mm, kk) == 1
                                p_conf_given_m = squeeze(p1_conf(mm,:,:));
                            elseif confResp_A_V(mm, kk) == 2
                                p_conf_given_m = squeeze(p2_conf(mm,:,:));
                            elseif confResp_A_V(mm, kk) == 3
                                p_conf_given_m = squeeze(p3_conf(mm,:,:));
                            elseif confResp_A_V(mm, kk) == 4
                                p_conf_given_m = squeeze(p4_conf(mm,:,:));
                            end
                            nLL_bimodal = nLL_bimodal - log(sum(sum(p_conf_given_m.*p_mAmV_given_sAsV)));
                        end
                    end
                elseif strcmp(model.strategy_loc,'MS') %Model selection
                    for mm = 1:num_modality
                        for kk = 1:num_rep

                            % confidence probability
                            if confResp_A_V(mm, kk) == 1
                                p_conf_given_m = squeeze(p1_conf(mm,:,:));
                            elseif confResp_A_V(mm, kk) == 2
                                p_conf_given_m = squeeze(p2_conf(mm,:,:));
                            elseif confResp_A_V(mm, kk) == 3
                                p_conf_given_m = squeeze(p3_conf(mm,:,:));
                            elseif confResp_A_V(mm, kk) == 4
                                p_conf_given_m = squeeze(p4_conf(mm,:,:));
                            end
                            nLL_bimodal = nLL_bimodal - log(sum(p_conf_given_m.*p_mAmV_given_sAsV));
                        end
                    end

                end
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

end