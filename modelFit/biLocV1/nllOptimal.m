function out = nllOptimal(sigA, sigV1, sigV2, sigP, pC1, c, model, data)

% fix parameters
muP = 0;
lapse = 0.01;

switch model.mode

    case 'optimize'

        f_logit = @(x) 1 ./ (1 + exp(-x));

        % combination of auditory and visual locations
        transA = model.sA * data.coefsA(2) + data.coefsA(1);
        transV = transA;
        sAV  = combvec(transA, transV);
        
        % fixed parameters
        num_a_stim = numel(transA);
        num_stim = size(sAV, 2);
        num_rep = model.num_rep;
        sigVs = [sigV1, sigV2];
        sigmaM = data.sigM;
        x = model.x;

        [loc_LL, conf_LL] = deal(NaN(1,numel(sigVs)));

        % loop by visual variable
        for i = 1:numel(sigVs)
            
            sigV = sigVs(i);

            % auditory locations (4) x visual locations (4) x postcues (2)
            % x rep
            data_resp = squeeze(data.org_resp(:,:,:,i,:));
            data_conf = squeeze(data.org_conf(:,:,:,i,:));

            mA                          = randn(num_stim, num_rep).*sigA + repmat(sAV(1,:)',[1, num_rep]);
            mV                          = randn(num_stim, num_rep).*sigV + repmat(sAV(2,:)',[1, num_rep]);

            mA                          = bounded(mA,min(x),max(x));
            mV                          = bounded(mV,min(x),max(x));

            JA                          = sigA^2;
            JP                          = sigP^2;

            JV                          = sigV^2;
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
            post_C1                     = pC1.*L_C1./(pC1.*L_C1 + (1-pC1).*L_C2);
            %posterior of separate causes = 1 - post_C1
            post_C2                     = 1 - post_C1;

            %compute the two intermediate location estimates
            %An integrated intermediate estimate is the sum of mA, mV and muP with
            %each weighted by their relative reliabilities
            %Eq. 12 in Körding et al., 2007
            sHat_C1                     = (mA./JA + mV./JV + muP/JP)./(1/JV + 1/JA + 1/JP);

            % shat_C1 = bounded(shat_C1,min(x),max(x));
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
            sHat_A = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
            sHat_V = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;
            loc(:,:,1,:) = reshape(sHat_A, [num_a_stim, num_a_stim, num_rep]);
            loc(:,:,2,:) = reshape(sHat_V, [num_a_stim, num_a_stim, num_rep]);

            % this model bases confidence on the variance of the full
            % posterior
            var_A = post_C1'.* 1/(1/JA + 1/JP) + post_C2'.* 1/(1/JV + 1/JA + 1/JP) + post_C1'.* post_C2' .* (sHat_A_C2' - sHat_C1').^2;
            var_V = post_C1'.* 1/(1/JV + 1/JP) + post_C2'.* 1/(1/JV + 1/JA + 1/JP) + post_C1'.* post_C2' .* (sHat_V_C2' - sHat_C1').^2;
            var(:,:,1,:) = reshape(var_A, [num_a_stim, num_a_stim, num_rep]);
            var(:,:,2,:) = reshape(var_V, [num_a_stim, num_a_stim, num_rep]);

            % convert variance to confidence variable
            conf = f_logit(var);
            pred_conf = conf > c;

            % likelihood
            p_loc = norm_dst(data_resp, loc, sigmaM, 1e-20);
            loc_LL(i) = sum(log(p_loc),'all');
            conf_LL(i) = sum(log( xor(pred_conf,data_conf) .* lapse + ~xor(pred_conf,data_conf) .* (1-lapse) ) ,'all');

        end

        out = - nansum(loc_LL(:)) - nansum(conf_LL(:));

    case 'predict'


end

    function x = bounded(x,LB,UB)
        x = max(LB,min(x, UB));
    end

    function p = norm_dst(x,mu,sigma,t)
        p = 1/sqrt(2*pi*sigma^2).*exp(-(x-mu).^2./(2*sigma^2)) + t;
    end

end