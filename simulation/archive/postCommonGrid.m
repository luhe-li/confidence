function posteriorC1 = postCommonGrid(sigA, sigV, muP, sigP,pcommon,x)
    % Range of the stimulus and measurement space
    range = x;
    
    % Pre-allocate probability matrices for C=1 and C=2
    likelihoodC1 = zeros(length(range), length(range));
    likelihoodC2 = zeros(length(range), length(range));

    JA     = sigA^2;
    JV     = sigV^2;
    JP     = sigP^2;
    const1 = JA*JV + JA*JP + JV*JP;
    constA = JA + JP;
    constV = JV + JP;
    % Calculate the likelihoods for C=1 and C=2 over the discrete space
    for i = 1:length(range)
        for j = 1:length(range)
            mV = range(i);
            mA = range(j); 
            
            %Eq. 4 in Körding et al., 2007
            likelihoodC1(i, j) = 1/(2*pi*sqrt(const1)).*exp(-0.5.*((mA - mV).^2 .* ...
            JP + (mV - muP).^2.*JA + (mA - muP).^2.* JV)./const1); 
            %Eq. 6 in Körding et al., 2007
            likelihoodC2(i, j) = 1/(2*pi*sqrt(constA*constV)).*exp(-0.5.*(mA - muP).^2./constA +...
            (mV - muP).^2./constV); 
        end
    end
    
    % Normalize the likelihoods so that they sum to 1
    likelihoodC1 = likelihoodC1 / sum(likelihoodC1, 'all');
    likelihoodC2 = likelihoodC2 / sum(likelihoodC2, 'all');
    
    % Calculate the posterior probability for C=1
    posteriorC1 = (likelihoodC1 * pcommon) ./ ...
                  ((likelihoodC1 * pcommon) + (likelihoodC2 * (1 - pcommon)));
end
