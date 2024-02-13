%This function simulates localization responses when a bimodal audiovisual
%stimulus sAV is presented. 

%input arguments:
%pCommon   : the common-cause prior about vision and audition
%nT        : number of simulated trials
%sA        : the location of the auditory component of the AV pair
%sV        : the location of the visual component of the AV pair
%aA        : the parameter that controls location-dependent biases in
%               auditory spatial perception
%bA        : the parameter that controls location-independent biases in
%               auditory spatial perception
%sigA      : measurement noise given an auditory stimulus
%sigV      : measurement noise given a visual stimulus
%muP       : the center of the prior distribution over stimulus location
%sigP      : the width of the prior distribution over stimulus location

%output arguments:
%r         : simulated localization responses (a matrix with size = 
%               3 decision strategies x 2 modalities x nT trials)
function [r,pdf] = simResp_CI_solutions_Multi(pCommon,nT, sA, sV, aA, bA, sigA, sigV, muP, sigP,x)
    %----------------------------------------------------------------------
    % YOUR CODE STARTS HERE
    %first initialize r
    %3 decision strategies (1. averaging, 2. selection, 3. matching)
    %2 modalities (1. auditory, 2. visual)
    %nT trials
    r     = NaN(3,2,nT); 

    %simulate measurements, which are drawn from Gaussian distributions
    % stochasticity starts here
    mA    = randn(1, nT).*sigA + (sA * aA + bA); 
    mV    = randn(1, nT).*sigV + sV;  
    
    %compute constants (these will come in handy when writing the equations
    %for the likelihood/posterior of a common cause and separate causes.
    JA     = sigA^2;
    JV     = sigV^2;
    JP     = sigP^2;
    const1 = JA*JV + JA*JP + JV*JP;
    constA = JA + JP;
    constV = JV + JP;

    %calculate the likelihood of a common cause and separate causes
    %Eq. 4 in Körding et al., 2007
    L_C1 = 1/(2*pi*sqrt(const1)).*exp(-0.5.*((mA - mV).^2 .* ...
            JP + (mV - muP).^2.*JA + (mA - muP).^2.* JV)./const1); 
    %Eq. 6 in Körding et al., 2007
    L_C2 = 1/(2*pi*sqrt(constA*constV)).*exp(-0.5.*(mA - muP).^2./constA +...
            (mV - muP).^2./constV); 

    %calculate posterior of a common cause and separate causes
    %Eq. 2 in Körding et al., 2007
    post_C1 = pCommon.*L_C1./(pCommon.*L_C1 + (1-pCommon).*L_C2); 
    %posterior of separate causes = 1 - post_C1
    post_C2 = 1 - post_C1;

    %compute the two intermediate location estimates
    %An integrated intermediate estimate is the sum of mA, mV and muP with
    %each weighted by their relative reliabilities
    %Eq. 12 in Körding et al., 2007
    shat_C1   = (mA./JA + mV./JV + muP/JP)./(1/JV + 1/JA + 1/JP); 
    %A segregated intermediate estimate is the sum of mA/mV and muP with
    %each weighted by their relative reliabilities
    %Eq. 11 in Körding et al., 2007
    sHat_A_C2 = (mA./JA + muP/JP)./(1/JA + 1/JP);
    sHat_V_C2 = (mV./JV + muP/JP)./(1/JV + 1/JP); 

    pdf.sHat_C1 = normpdf(repmat(x,nT,1),shat_C1',1/sqrt(1/JV + 1/JA + 1/JP));
    pdf.sHat_A_C2 = normpdf(repmat(x,nT,1),sHat_A_C2',1/sqrt(1/JA + 1/JP));
    pdf.sHat_V_C2 = normpdf(repmat(x,nT,1),sHat_V_C2',1/sqrt(1/JV + 1/JP));
    pdf.AMatching = pdf.sHat_C1 .* post_C1' + pdf.sHat_A_C2 .* post_C2';
    pdf.VMatching = pdf.sHat_C1 .* post_C1' + pdf.sHat_V_C2 .* post_C2';
    pdf.Aaveraging = normpdf(repmat(x,nT,1),(post_C1.* shat_C1 + post_C2.* sHat_A_C2)', ...
        sqrt( post_C1.* 1/(1/JV + 1/JA + 1/JP) + post_C2.* 1/(1/JA + 1/JP) + post_C1 .* post_C2 .* (shat_C1 - sHat_A_C2).^2)');
    pdf.Vaveraging = normpdf(repmat(x,nT,1),(post_C1.* shat_C1 + post_C2.* sHat_V_C2)', ...
        sqrt( post_C1.* 1/(1/JV + 1/JA + 1/JP) + post_C2.* 1/(1/JV + 1/JP) + post_C1 .* post_C2 .* (shat_C1 - sHat_V_C2).^2)');
    % Matching is effectively conditional probability
    
    % standardize the pdfs
    pdf.sHat_C1 = pdf.sHat_C1 ./ sum(pdf.sHat_C1,2);
    pdf.sHat_A_C2 = pdf.sHat_A_C2 ./ sum(pdf.sHat_A_C2,2);
    pdf.sHat_V_C2 = pdf.sHat_V_C2 ./ sum(pdf.sHat_V_C2,2);
    pdf.AMatching = pdf.AMatching ./ sum(pdf.AMatching,2);
    pdf.VMatching = pdf.VMatching ./ sum(pdf.VMatching,2);
    pdf.Aaveraging = pdf.Aaveraging ./ sum(pdf.Aaveraging,2);
    pdf.Vaveraging = pdf.Vaveraging ./ sum(pdf.Vaveraging,2);
    %compute the final location estimates if we assume model averaging.
    %Based on this strategy, the final location estimate is the sum of the
    %two intermediate location estimates, weighted by the corresponding
    %causal structure.
    %Eq. 4 in Wozny et al., 2010
    r(1,1,:) = post_C1.* shat_C1 + post_C2.* sHat_A_C2;
    r(1,2,:) = post_C1.* shat_C1 + post_C2.* sHat_V_C2;
    
    %compute the final location estimates if we assume model selection.
    %Based on this strategy, the final location estimate depends purely on
    %the causal structure that is more probable.
    %Eq. 5 in Wozny et al., 2010
    r(2,1:2,:)         = repmat(shat_C1,[2 1]);
    r(2,1,post_C1<0.5) = sHat_A_C2(post_C1<0.5);
    r(2,2,post_C1<0.5) = sHat_V_C2(post_C1<0.5);
    pdf.selectionComm = ~(post_C1<0.5);

    %compute the final location estimates if we assume probability matching.
    %Based on this strategy, the final location estimate is the integrated
    %one with a probability of post_C1, and is the segregated one with a
    %probability of post_C2.
    %Eq. 6 in Wozny et al., 2010
    idx           = rand(1,nT);
    idx_C2        = (idx > post_C1);
    r(3,1:2,:)    = repmat(shat_C1,[2 1]);
    r(3,1,idx_C2) = sHat_A_C2(idx_C2);
    r(3,2,idx_C2) = sHat_V_C2(idx_C2);
    pdf.matchingComm = ~idx_C2;
    %----------------------------------------------------------------------
end


