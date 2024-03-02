
clear; clc; close all;

%% set up

sub = 1;

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(project_dir,'exptCode','biloc');
addpath(genpath(fullfile(project_dir, 'data','uniloc')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% organize data

% load V first to avoid overwriting sortedResp
load(sprintf('uniLoc_sub%i_ses%s', sub, '-V'))
load(sprintf('uniLoc_sub%i_ses%s', sub, '-A'))
% set up
nRep = ExpInfo.nRep;
targIdx = unique([ExpInfo.randAudIdx]);
targDegs = unique(ExpInfo.speakerLocVA(ExpInfo.randAudIdx));
targPx = unique(ExpInfo.speakerLocPixel(ExpInfo.randAudIdx));
targNum = length(targIdx);
rA = NaN(targNum,nRep);
rV1 = NaN(targNum,nRep);
rV2 = NaN(targNum,nRep);

for i = 1:targNum
    targInd   = targIdx(i);
    estDegsA  = [sortedResp.response_deg];
    estDegsV1 = [sortedReli1Resp.response_deg];
    estDegsV2 = [sortedReli2Resp.response_deg];

    rA(i,:) = estDegsA([sortedResp.target_idx] == targInd);
    rV1(i,:) = estDegsV1([sortedReli1Resp.target_idx] == targInd);
    rV2(i,:) = estDegsV2([sortedReli2Resp.target_idx] == targInd);
end

%% fit a linear line to both

x = repmat(targDegs',1,nRep);
x2 = repmat(targDegs',1,nRep*2);

mdlA = fitlm(x(:),rA(:));
coefsA = table2array(mdlA.Coefficients(:,1));
fitRA = x(:,1) .* coefsA(2) + coefsA(1);

mdlV = fitlm(x2(:),[rV1(:);rV2(:)]);
coefsV = table2array(mdlV.Coefficients(:,1));
fitSV = (fitRA - coefsV(1)) ./ coefsV(2);

Transfer.degCoeff(1, :) = coefsA;
Transfer.degCoeff(2, :) = coefsV;

%%
rA = NaN(targNum,nRep);
rV1 = NaN(targNum,nRep);
rV2 = NaN(targNum,nRep);

for i = 1:targNum
    targInd   = targIdx(i);
    estPxA  = [sortedResp.response_pixel];
    estPxV1 = [sortedReli1Resp.response_pixel];
    estPxV2 = [sortedReli2Resp.response_pixel];

    rA(i,:) = estPxA([sortedResp.target_idx] == targInd);
    rV1(i,:) = estPxV1([sortedReli1Resp.target_idx] == targInd);
    rV2(i,:) = estPxV2([sortedReli2Resp.target_idx] == targInd);
end

x = repmat(targPx',1,nRep);
x2 = repmat(targPx',1,nRep*2);

mdlA = fitlm(x(:),rA(:));
coefsA = table2array(mdlA.Coefficients(:,1));
fitRA = x(:,1) .* coefsA(2) + coefsA(1);

mdlV = fitlm(x2(:),[rV1(:);rV2(:)]);
coefsV = table2array(mdlV.Coefficients(:,1));
fitSV = (fitRA - coefsV(1)) ./ coefsV(2);

Transfer.PxCoeff(1, :) = coefsA;
Transfer.PxCoeff(2, :) = coefsV;

save(fullfile(out_dir, sprintf('AVbias_sub%i',sub)), 'Transfer')

% %%%%%%%%%%%%%%%%%%%%%%%%%% check plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% plot response to check
% 
% figure
% set(gcf, 'Position',[0 0 900 300])
% 
% sgtitle(sprintf('S%i',sub_add))
% 
% subplot(1,3,1)
% scatter(targDegs,rA,'filled','r')
% hold on
% plot(-30:30,-30:30,'--')
% hold off
% xlim([-30,30])
% xlabel('Stimulus')
% ylabel('Estimation')
% title('A')
% 
% 
% subplot(1,3,2)
% scatter(targDegs,rV1,'filled','r')
% hold on
% plot(-30:30,-30:30,'--')
% hold off
% xlim([-30,30])
% xlabel('Stimulus')
% title('V1')
% 
% subplot(1,3,3)
% scatter(targDegs,rV2,'filled','r')
% hold on
% plot(targDegs, fitSV2,'k-')
% plot(-30:30,-30:30,'--')
% hold off
% xlim([-30,30])
% xlabel('Stimulus')
% title('V2')

% %% optinal: check visual response, V2 vs V1
% figure
% axis equal
% axis square
% scatter(rA(:),rV2(:),'filled','r')
% xlim([-35,35])
% ylim([-35,35])
