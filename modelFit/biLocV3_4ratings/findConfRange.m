function [lb, ub] = findConfRange(aA, bA, sigA, sigV1, sigV2, sigP, pCommon, fixP)
% findConfRange takes in free parameters other than criteria and calculates
% the plausible range for variance/uncertainty.

% 24/05 LL

sigVs = [sigV1, sigV2];
num_rep = fixP.num_rep;
muP = 0;

%% simulate variance for all location combinations and visual reliabilities
for aa = 1:numel(fixP.sA)

    sA = fixP.sA(aa);

    for  vv = 1:numel(fixP.sV)

        sV = fixP.sV(vv);

        for rr = 1:2

            sigV = sigVs(rr);

            if fixP.model_ind == 1

                variance(aa, vv, 1, rr, :) = repmat(sigA, [1, num_rep]); % normalized by sigA
                variance(aa, vv, 2, rr,: ) = repmat(sigV, [1, num_rep]); % normalized by sigV

            else

                %simulate measurements, which are drawn from Gaussian distributions
                % stochasticity starts here
                mA    = randn(1, num_rep).*sigA + (sA * aA + bA);
                mV    = randn(1, num_rep).*sigV + sV;

                %compute constants (these will come in handy when writing the equations
                %for the likelihood/posterior of a common cause and separate causes.
                JA    = sigA^2;
                JV    = sigV^2;
                JP    = sigP^2;
                const1= JA*JV + JA*JP + JV*JP;
                constA= JA + JP;
                constV= JV + JP;

                %calculate the likelihood of a common cause and separate causes
                %Eq. 4 in Körding et al., 2007
                L_C1  = 1/(2*pi*sqrt(const1)).*exp(-0.5.*((mA - mV).^2 .* ...
                    JP + (mV - muP).^2.*JA + (mA - muP).^2.* JV)./const1);
                %Eq. 6 in Körding et al., 2007
                L_C2  = 1/(2*pi*sqrt(constA*constV)).*exp(-0.5.*(mA - muP).^2./constA +...
                    (mV - muP).^2./constV);

                %calculate posterior of a common cause and separate causes
                %Eq. 2 in Körding et al., 2007
                post_C1    = pCommon.*L_C1./(pCommon.*L_C1 + (1-pCommon).*L_C2);
                %posterior of separate causes
                post_C2    = 1 - post_C1;

                if fixP.model_ind == 2

                    JA = sigA^2;
                    JV = sigV^2;
                    JP = sigP^2;
                    variance(aa, vv, 1:2, rr, :)= repmat(1/(1/JV + 1/JA + 1/JP), [2, num_rep]);
                    variance(aa, vv, 1, rr, post_C1<0.5)  = 1/(1/JA + 1/JP)./sigA;
                    variance(aa, vv, 2, rr, post_C1<0.5)  = 1/(1/JV + 1/JP)./sigV;

                else
                    %compute the two intermediate location estimates
                    %An integrated intermediate estimate is the sum of mA, mV and muP with
                    %each weighted by their relative reliabilities
                    %Eq. 12 in Körding et al., 2007
                    sHat_C1    = (mA./JA + mV./JV + muP/JP)./(1/JV + 1/JA + 1/JP);

                    %A segregated intermediate estimate is the sum of mA/mV and muP with
                    %each weighted by their relative reliabilities
                    %Eq. 11 in Körding et al., 2007
                    sHat_A_C2  = (mA./JA + muP/JP)./(1/JA + 1/JP);
                    sHat_V_C2  = (mV./JV + muP/JP)./(1/JV + 1/JP);

                    %compute the final location estimates if we assume model averaging.
                    %Based on this strategy, the final location estimate is the sum of the
                    %two intermediate location estimates, weighted by the corresponding
                    %causal structure.
                    %Eq. 4 in Wozny et al., 2010
                    shat(1,:) = post_C1.* sHat_C1 + post_C2.* sHat_A_C2;
                    shat(2,:) = post_C1.* sHat_C1 + post_C2.* sHat_V_C2;

                    variance(aa, vv, 1, rr, :)= (post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JA + 1/JP) + post_C1.* post_C2 .* (sHat_A_C2 - sHat_C1).^2)./sigA;
                    variance(aa, vv, 2, rr, :)= (post_C1./(1/JV + 1/JA + 1/JP) + post_C2./(1/JV + 1/JP) + post_C1.* post_C2 .* (sHat_V_C2 - sHat_C1).^2)./sigV;

                end
            end
        end

    end

end

%% oganize by discrepancy

% {diff} cue x reliability x rep
[var_by_diff, all_diffs] = org_by_diffs(variance, fixP.sA);

% % visualize variance by discrepancy
% figure; hold on
% tiledlayout(2,2);
all_var = [];
for cue = 1:2
    for rel = 1: numel(sigVs)
        %         nexttile
        %         hold on
        %         xlim([0,3])
        for diff = 1:numel(all_diffs)
            i_var = squeeze(var_by_diff{diff}(cue, rel, :))';
            all_var = [all_var, i_var];
            %             histogram(i_var,'BinWidth',0.05)
        end

    end
    %     legend
end

lb = min(all_var,[],'all');
ub = max(all_var,[],'all');

end