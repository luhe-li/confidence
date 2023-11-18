function [d_reshaped, d_reshaped_mean, d_reshaped_SD, numT_AV] = reshapeResp(d,nT_perPair)
%reshape the data from lenS x lenS to lenS x lenD
%-------------------------Spatial discrepancy(V-A)-----------------------------
%       -24,       -16,       -8,        0,         8,        16,        24
%------------------------------------------------------------------------------
%V1: p(_|A4,V1) p(_|A3,V1)  p(_|A2,V1)  p(_|A1,V1)     NaN       NaN       NaN
%V2:    NaN     p(_|A4,V2)  p(_|A3,V2)  p(_|A2,V2)  p(_|A1,V2)   NaN       NaN  
%V3:    NaN         NaN     p(_|A4,V3)  p(_|A3,V3)  p(_|A2,V3) p(_|A1,V3)  NaN   
%V4:    NaN         NaN        NaN      p(_|A4,V4)  p(_|A3,V4) p(_|A2,V4) p(_|A1,V4)
global lenC lenP lenS lenD lenM
d_reshaped = NaN(lenC, lenP, lenS, lenD);
[d_reshaped_mean, d_reshaped_SD] = deal(NaN(lenC, lenP, lenD));
%mean p(reporting 'C=1') as a function of spatial discrepancy
for i = 1:lenC
    for j = 1:lenP
        for l = 1:lenS %for each visual/auditory stimulus location
            idx_lb = l; idx_ub = floor(lenD/2)+l;
            %need to use fliplr.m because the order is 
            %p('C=1'|A4,V1), p('C=1'|A3,V1), p('C=1'|A2,V1), p('C=1'|A1,V1)
            d_reshaped(i,j,l,idx_lb:idx_ub) = fliplr(squeeze(d(i,j,:,l))');
        end
        %calculate the mean percentage of reporting C=1 as a function of discrepancy
        d_reshaped_mean(i,j,:) = squeeze(nanmean(d_reshaped(i,j,:,:)));
        if i==1 && j==1 
            numT_AV = nT_perPair*lenM.*sum(~isnan(squeeze(d_reshaped(i,j,:,:)))); 
        end
        d_reshaped_SD(i,j,:) = calcSD(squeeze(d_reshaped_mean(i,j,:))', numT_AV); 
    end
end

function sd = calcSD(p,n)%define a function that calculates the SD for bernoulli distribution
sd = sqrt(p.*(1-p)./n);

