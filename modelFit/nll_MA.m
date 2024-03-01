function out = nll_MA(aA, bA, sig_AV_A, sig_AV_V1, sig_AV_V2 sigP, pC1, c, data, model)

mu_P = 0;
aV = 1;
bV = 0;

switch mode
    
    case 'init'



    case 'fit'

        %% bimodal
        x = fixP.x;
        num_rep = fixP.num_rep;
        loc = NaN(2,num_rep);
        p_conf = NaN(2, num_rep);

        mA                          = randn(1, num_rep).*sigA + (sA * aA + bA);
        mV                          = randn(1, num_rep).*sigV + sV;

        mA                          = bounded(mA,min(x),max(x));
        mV                          = bounded(mV,min(x),max(x));

        JA                          = sigA^2;
        JV                          = sigV^2;
        JP                          = sigP^2;
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
        post_C1                     = pCommon.*L_C1./(pCommon.*L_C1 + (1-pCommon).*L_C2);
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
        loc(1,:) = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
        loc(2,:) = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;



        var(1,:) = post_C1'.* 1/(1/JA + 1/JP) + post_C2'.* 1/(1/JV + 1/JA + 1/JP) + post_C1'.* post_C2' .* (sHat_A_C2' - sHat_C1').^2;
        var(2,:) = post_C1'.* 1/(1/JV + 1/JP) + post_C2'.* 1/(1/JV + 1/JA + 1/JP) + post_C1'.* post_C2' .* (sHat_V_C2' - sHat_C1').^2;

        p_conf = var > c;

        LL = nt_yes * log()' + nt_no * log()';
        out = -LL;

    case 'predict'


end


end