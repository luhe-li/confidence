function sim_binaryJdg = sim_binaryJdg(sim_pC1, nT_perPair)

global lenC lenP lenS
sim_binaryJdg = NaN(lenC, lenP, lenS, lenS, nT_perPair); %boolean
for i = 1:lenC
    for j = 1:lenP
        for l = 1:lenS
            for m = 1:lenS
                sim_binaryJdg(i,j,l,m,:) = getBinary(sim_pC1(i,j,l,m), nT_perPair);
            end
        end
    end
end

function binaryResp = getBinary(p,n)
binaryResp = [ones(1,round(p*n)), zeros(1,round((1-p)*n))];

