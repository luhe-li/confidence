function out = nllOptimal(sigA, sigV1, sigV2, sigP, pC1, c, data, model)

% fix parameters
aA = 1;
bA = 0;
muP = 0;
lapse = 0.01;

switch mode

    case 'fit'

        x = model.x;
        num_rep = model.num_rep;
        num_stim = numel(sA);
        sigVs = [sigV1, sigV2];
        transA = model.sA * aA + bA;
        transV = transA;

        for i = 1:numel(sigVs)

            sigV = sigVs(i);
            mA                          = randn(num_rep, num_stim).*sigA + transA;
            mV                          = randn(num_rep, num_stim).*sigV + transV;

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
            loc(1,:,:) = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
            loc(2,:,:) = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;

            % base confidence on variance
            var(1,:,:) = post_C1'.* 1/(1/JA + 1/JP) + post_C2'.* 1/(1/JV + 1/JA + 1/JP) + post_C1'.* post_C2' .* (sHat_A_C2' - sHat_C1').^2;
            var(2,:,:) = post_C1'.* 1/(1/JV + 1/JP) + post_C2'.* 1/(1/JV + 1/JA + 1/JP) + post_C1'.* post_C2' .* (sHat_V_C2' - sHat_C1').^2;
            p_conf = var > c;
            p_conf_lapsed = NaN(size(p_conf));
            p_conf_lapsed(p_conf == 1) = 1 - lapse;
            p_conf_lapsed(p_conf == 0) = lapse;

            % likelihood
            p_loc = norm_dst(loc_resp, loc, sigmaM, 1e-20);
            loc_LL = log(p_loc);
            conf_LL = nt_yes * log()' + nt_no * log()';

        end

        out = -conf_LL;

    case 'predict'


end

    function x = bounded(x,LB,UB)
        x = max(LB,min(x, UB));
    end

    function p = norm_dst(x,mu,sigma,t)
        p = 1/sqrt(2*pi*sigma^2).*exp(-(x-mu).^2./(2*sigma^2)) + t;
    end

end