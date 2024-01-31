
stds = NaN(3,8);
targInds = unique([sortedResp.target_idx]);
targDegs = unique([sortedResp.target_deg]);
estA = NaN(8,20);
estV1 = NaN(8,20);
estV2 = NaN(8,20);

for i = 1:8
    targInd = targInds(i);
    estDegsA = [sortedResp.response_deg];
    estDegsV1 = [sortedReli1Resp.response_deg];
    estDegsV2 = [sortedReli2Resp.response_deg];

    estA(i,:) = estDegsA([sortedResp.target_idx] == targInd);
    estV1(i,:) = estDegsV1([sortedReli1Resp.target_idx] == targInd);
    estV2(i,:) = estDegsV1([sortedReli2Resp.target_idx] == targInd);

    
end

%%
figure
subplot(1,3,1)
scatter(targDegs,estA,'filled','r')
hold on
plot(-30:30,-30:30,'--')
hold off
xlim([-30,30])
xlabel('Stimulus')
ylabel('Estimation')
title('A')


subplot(1,3,2)
scatter(targDegs,estV1,'filled','r')
hold on
plot(-30:30,-30:30,'--')
hold off
xlim([-30,30])
xlabel('Stimulus')
title('V1')

subplot(1,3,3)
scatter(targDegs,estV2,'filled','r')
hold on
plot(-30:30,-30:30,'--')
hold off
xlim([-30,30])
xlabel('Stimulus')
title('V2')


%%




