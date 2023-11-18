function [r, counts] = sim_locResp(prob_r, s, numT)
global lenC lenP lenM lenS
r = zeros(lenC, lenP, lenS, lenS, lenM, numT);
counts = zeros(lenC, lenP, lenS, lenS, lenM, length(s));
for i = 1:lenC
    for j = 1:lenP
        for k = 1:lenS
            for l = 1:lenS
                for m = 1:lenM
                    prob_r_temp = squeeze(prob_r(i,j,k,l,m,:))';
                    prob_r_norm = prob_r_temp./sum(prob_r_temp);
                    prob_r_cum = [0,cumsum(prob_r_norm)];
                    randN = rand(1,numT);
                    bool_less = (randN' < prob_r_cum);
                    for n = 1:numT
                        idx = find(bool_less(n,:)==1, 1);
                        counts(i,j,k,l,m,idx) = counts(i,j,k,l,m,idx) + 1;
                        r(i,j,k,l,m,n) = s(idx);
                    end
                end
            end
        end
    end
end

end