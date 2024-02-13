
nRep = ExpInfo.nRep;
targInds = unique([sortedResp.target_idx]);
targDegs = unique([sortedResp.target_deg]);
targNum = length(targInds);
rA = NaN(targNum,nRep);
rV1 = NaN(targNum,nRep);
rV2 = NaN(targNum,nRep);

for i = 1:targNum
    targInd = targInds(i);
    estDegsA = [sortedResp.response_deg];
    estDegsV1 = [sortedReli1Resp.response_deg];
    estDegsV2 = [sortedReli2Resp.response_deg];

    rA(i,:) = estDegsA([sortedResp.target_idx] == targInd);
    rV1(i,:) = estDegsV1([sortedReli1Resp.target_idx] == targInd);
    rV2(i,:) = estDegsV2([sortedReli2Resp.target_idx] == targInd);
end

%%
figure
subplot(1,3,1)
scatter(targDegs,rA,'filled','r')
hold on
plot(-30:30,-30:30,'--')
hold off
xlim([-30,30])
xlabel('Stimulus')
ylabel('Estimation')
title('A')


subplot(1,3,2)
scatter(targDegs,rV1,'filled','r')
hold on
plot(-30:30,-30:30,'--')
hold off
xlim([-30,30])
xlabel('Stimulus')
title('V1')

subplot(1,3,3)
scatter(targDegs,rV2,'filled','r')
hold on
plot(-30:30,-30:30,'--')
hold off
xlim([-30,30])
xlabel('Stimulus')
title('V2')

%%
scatter(rA(:),rV2(:),'filled','r')
xlim([-35,35])
ylim([-35,35])

%%
x = repmat(targDegs',1,nRep);

mdlA = fitlm(x(:),rA(:));
coefsA = table2array(mdlA.Coefficients(:,1));
fitRA = x(:,1) .* coefsA(2) + coefsA(1);

mdlV1 = fitlm(x(:),rV1(:));
coefsV1 = table2array(mdlV1.Coefficients(:,1));
fitSV1 = (fitRA - coefsV1(1)) ./ coefsV1(2);

mdlV2 = fitlm(x(:),rV2(:));
coefsV2 = table2array(mdlV2.Coefficients(:,1));
fitSV2 = (fitRA - coefsV2(1)) ./ coefsV2(2);



