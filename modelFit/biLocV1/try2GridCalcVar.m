clear all
speaker_span       = 65.5 * 2; % cm
sitting_dist       = 113; % cm
screen_width       = 170; % cm
screenX            = 1024; % pixel
screenXdeg         = rad2deg(atan(screen_width / 2 / sitting_dist)) .* 2;
screen_mid         = screenX ./2;
num_speaker_int    = 15; % 15 intervals between 16 speakers
cm_per_aud_ind     = speaker_span / num_speaker_int;
pixel_per_cm       = screenX / screen_width;
aud_level          = [6 8 9 11];
aud_VA             = -30:4:30;
deg_per_px         = screenXdeg / screenX;

screenX       = screenXdeg;
x             = -screenXdeg /2 : deg_per_px : screenXdeg/2;

%%
sigA = 3;
sigP = 20;
sigV = 5;

JA                          = sigA^2;
JP                          = sigP^2;

JV                          = sigV^2;
const1                      = JA*JV + JA*JP + JV*JP;
constA                      = JA + JP;
constV                      = JV + JP;
pC1 = 0.57;
muP = 0;
sA = 2;
sV = 2;
%%
for aInd = 1:length(x)
    for vInd = 1:length(x)
        mA = x(aInd);
        mV = x(vInd);


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


        % this model bases confidence on the variance of the full
        % posterior
        var_A(aInd,vInd) = post_C1'.* 1/(1/JA + 1/JP) + post_C2'.* 1/(1/JV + 1/JA + 1/JP) + post_C1'.* post_C2' .* (sHat_A_C2' - sHat_C1').^2;
        var_V(aInd,vInd) = post_C1'.* 1/(1/JV + 1/JP) + post_C2'.* 1/(1/JV + 1/JA + 1/JP) + post_C1'.* post_C2' .* (sHat_V_C2' - sHat_C1').^2;
        pdfAV(aInd,vInd) = normpdf(mA,sA,sigA) * normpdf(mV,sV,sigV);
    end
end

%%
pdfAV = pdfAV ./ sum(pdfAV,"all");
pdfAV = pdfAV(:);
%%
[N,edges,bin] = histcounts(var_V);
binnedPDF = NaN(1,length(N));

for i = 1:length(N)
    binnedPDF(i) = sum(pdfAV(bin == i));
end
plot(edges(1:end-1),N .* binnedPDF,'-o')







