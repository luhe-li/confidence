subjNs = [3,4,5,6,8,9,11,12,13,15,16,17,18,19,20];
subjIs = {'PW','SW','HL','YZ','NH','ZZ','BB','ZY','MR','AD','SM','SX','ZL','RE','MM'};
order_dict = [2,1;1,2;2,1;1,2;2,1;1,2;2,1;2,1;1,2;2,1;1,2;1,2;2,1;1,2;2,1];
cond   = {'congruent', 'incongruent'};
nS     = length(subjNs);
prop_C1 = NaN(nS, length(cond));
for i = 1:nS
    for j = 1:length(cond)
        sesN = order_dict(i,j);
        addpath(genpath(['/Users/hff/Desktop/NYU/Project2/Experiment code/',...
            'Adaptation v2/Data/',subjIs{i}]));
        try
            C = load(['Adaptation_',cond{sesN}, '_sub',num2str(subjNs(i)),...
                '_session', num2str(j), '.mat'], 'Adaptation_data');
            D = C.Adaptation_data{end}.unity;
        catch
            C = load(['Adaptation_',cond{sesN}, '_sub',num2str(subjNs(i)),...
                '_session', num2str(j), '.mat'], 'Adaptation_incongruent_data');  
            D = C.Adaptation_incongruent_data{end}.unity;
        end
        prop_C1(i,sesN) = nansum(D==1)./(length(D).*0.2);
    end
end

%%
disp(mean(prop_C1))
disp(std(prop_C1)./sqrt(nS))


