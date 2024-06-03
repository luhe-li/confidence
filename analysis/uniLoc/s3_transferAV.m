
clear; clc; close all;

%% set up

sub = 12;

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(project_dir,'exptCode','biloc');
addpath(genpath(fullfile(project_dir, 'data','uniloc')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% organize data

load(sprintf('uniLoc_sub%i_ses%s', sub, '-A'))
% set up
nRep = ExpInfo.nRep;
targIdx = unique([ExpInfo.randAudIdx]);
targDegs = unique(ExpInfo.speakerLocVA(ExpInfo.randAudIdx));
targPx = unique(ExpInfo.speakerLocPixel(ExpInfo.randAudIdx));
targNum = length(targIdx)-2;
rA = NaN(targNum,nRep);
for i = 1:targNum
    targInd   = targIdx(i+1);
    estDegsA  = [sortedResp.response_deg];
    rA(i,:) = estDegsA([sortedResp.target_idx] == targInd);
end

%% fit a linear line to both

x = repmat(targDegs(2:targNum+1)',1,nRep);

mdlA = fitlm(x(:),rA(:));
coefsA = table2array(mdlA.Coefficients(:,1));

Transfer.degCoeff(1, :) = coefsA;

%%
rA = NaN(targNum,nRep);

for i = 1:targNum
    targInd   = targIdx(i+1);
    estPxA  = [sortedResp.response_pixel];

    rA(i,:) = estPxA([sortedResp.target_idx] == targInd);
end

x = repmat(targPx(2:targNum+1)',1,nRep);

mdlA = fitlm(x(:),rA(:));
coefsA = table2array(mdlA.Coefficients(:,1));
fitRA = x(:,1) .* coefsA(2) + coefsA(1);
fitSV = fitRA - 512;

Transfer.PxCoeff(1, :) = coefsA;

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
