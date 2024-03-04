function out = nllHeuristic(freeParam, model, data)

switch model.mode

    case 'initiate'

        out.paraID                   = {'aA','bA','\sigma_{A}','\sigma_{V1}','\sigma_{V2}','\sigma_{P}','p_{common}','criterion'};
        out.num_para                 = length(out.paraID);

        % hard bounds, the range for LB, UB, larger than soft bounds
        paraH.aA                     = [ 0.5,     2]; % degree
        paraH.bA                     = [  -8,     8]; % degree
        paraH.sigA                   = [1e-2,    10]; % degree
        paraH.sigV1                  = [1e-2,    10]; % degree
        paraH.sigV2                  = [1e-2,    10]; % degree
        paraH.sigP                   = [   1,    30]; % degree
        paraH.pC1                    = [1e-4,1-1e-4]; % weight
        paraH.c                      = [1e-4,1-1e-4]; % weight

        % soft bounds, the range for PLB, PUB
        paraS.aA                     = [ 0.8,   1.2]; % degree
        paraS.bA                     = [  -4,     4]; % degree
        paraS.sigA                   = [   4,     6]; % degree
        paraS.sigV1                  = [   2,     6]; % degree
        paraS.sigV2                  = [   4,     6]; % degree
        paraS.sigP                   = [   1,    20]; % degree
        paraS.pC1                    = [ 0.4,   0.6]; % weight
        paraS.c                      = [ 0.2,   0.5]; % weight

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
        out.init                     = getInit(out.lb, out.ub, model.num_runs, model.num_runs);

    case {'optimize','predict'}

        aA                           = freeParam(1);
        bA                           = freeParam(2);
        sigA                         = freeParam(3);
        sigV1                        = freeParam(4);
        sigV2                        = freeParam(5);
        sigP                         = freeParam(6);
        pC1                          = freeParam(7);
        c                            = freeParam(8);

        % logistic function that converts the confidence variable
        f_logit                      = @(x) 1 ./ (1 + exp(-x));

        % tentative parameters (can be fitted later)
        muP                          = 0;
        lapse                        = 0.01;

        % combination of auditory and visual locations
        transA                       = model.sA * aA + bA;
        transV                       = transA;
        sAV                          = combvec(transA, transV);

        % fixed parameters
        num_a_stim                   = numel(transA);
        num_stim                     = size(sAV, 2);
        num_rep                      = model.num_rep;
        sigVs                        = [sigV1, sigV2];
        sigmaM                       = data.sigM;
        x                            = model.x;

        [loc_LL, conf_LL]            = deal(NaN(1,numel(sigVs)));
        [all_loc, all_var, all_conf] = deal(NaN(num_a_stim,num_a_stim,2,2,num_rep));

        % loop by visual variable
        for i                        = 1:numel(sigVs)

            sigV                         = sigVs(i);

            % auditory locations (4) x visual locations (4) x postcues (2)
            % x rep
            data_resp                    = squeeze(data.org_resp(:,:,:,i,:));
            data_conf                    = squeeze(data.org_conf(:,:,:,i,:));

            mA                           = randn(num_stim, num_rep).*sigA + repmat(sAV(1,:)',[1, num_rep]);
            mV                           = randn(num_stim, num_rep).*sigV + repmat(sAV(2,:)',[1, num_rep]);

            mA                           = bounded(mA,min(x),max(x));
            mV                           = bounded(mV,min(x),max(x));

            JA                           = sigA^2;
            JP                           = sigP^2;
            JV                           = sigV^2;
            const1                       = JA*JV + JA*JP + JV*JP;
            constA                       = JA + JP;
            constV                       = JV + JP;

            %calculate the likelihood of a common cause and separate causes
            %Eq. 4 in Körding et al., 2007
            L_C1                         = 1/(2*pi*sqrt(const1)).*exp(-0.5.*((mA - mV).^2 .* ...
                JP + (mV - muP).^2.*JA + (mA - muP).^2.* JV)./const1);
            %Eq. 6 in Körding et al., 2007
            L_C2                         = 1/(2*pi*sqrt(constA*constV)).*exp(-0.5.*(mA - muP).^2./constA +...
                (mV - muP).^2./constV);

            %calculate posterior of a common cause and separate causes
            %Eq. 2 in Körding et al., 2007
            post_C1                      = pC1.*L_C1./(pC1.*L_C1 + (1-pC1).*L_C2);
            %posterior of separate causes = 1 - post_C1
            post_C2                      = 1 - post_C1;

            %compute the two intermediate location estimates
            %An integrated intermediate estimate is the sum of mA, mV and muP with
            %each weighted by their relative reliabilities
            %Eq. 12 in Körding et al., 2007
            sHat_C1                      = (mA./JA + mV./JV + muP/JP)./(1/JV + 1/JA + 1/JP);

            % shat_C1 = bounded(shat_C1,min(x),max(x));
            %A segregated intermediate estimate is the sum of mA/mV and muP with
            %each weighted by their relative reliabilities
            %Eq. 11 in Körding et al., 2007
            sHat_A_C2                    = (mA./JA + muP/JP)./(1/JA + 1/JP);
            sHat_V_C2                    = (mV./JV + muP/JP)./(1/JV + 1/JP);

            %compute the final location estimates if we assume model averaging.
            %Based on this strategy, the final location estimate is the sum of the
            %two intermediate location estimates, weighted by the corresponding
            %causal structure.
            %Eq. 4 in Wozny et al., 2010
            sHat_A                       = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
            sHat_V                       = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;
            loc(:,:,1,:)                 = reshape(sHat_A, [num_a_stim, num_a_stim, num_rep]);
            loc(:,:,2,:)                 = reshape(sHat_V, [num_a_stim, num_a_stim, num_rep]);

            % this model bases confidence on the variance of sensory noise
            var(:,:,1,:)                 = repmat(JA, [num_a_stim, num_a_stim, num_rep]);
            var(:,:,2,:)                 = repmat(JV, [num_a_stim, num_a_stim, num_rep]);

            % convert variance to confidence variable
            conf                         = f_logit(var);
            pred_conf                    = conf > c;

            % likelihood
            p_loc                        = norm_dst(data_resp, loc, sigmaM, 1e-20);
            loc_LL(i)                    = sum(log(p_loc),'all');
            conf_LL(i)                   = sum(log( xor(pred_conf,data_conf) .* lapse + ~xor(pred_conf,data_conf) .* (1-lapse) ) ,'all');

            all_loc(:,:,:,i,:)           = loc;
            all_var(:,:,:,i,:)           = var;
            all_conf(:,:,:,i,:)          = pred_conf;

        end

        if strcmp(model.mode, 'optimize')

            out                          = - nansum(loc_LL(:)) - nansum(conf_LL(:));

        elseif strcmp(model.mode, 'predict')

            out.resp                     = all_loc;
            out.var                      = all_var;
            out.conf                     = all_conf;

        end

end

    function x                   = bounded(x,LB,UB)
        x                            = max(LB,min(x, UB));
    end

    function p                   = norm_dst(x,mu,sigma,t)
        p                            = 1/sqrt(2*pi*sigma^2).*exp(-(x-mu).^2./(2*sigma^2)) + t;
    end

end