% v1: assume separate audiovisual bias for visual reliability 1 and 2

clear; clc; close all;

%% set up

sub = 1;

%% manage path

cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));
out_dir                          = fullfile(cur_dir, 's3AVbias');
addpath(genpath(fullfile(project_dir, 'data','uniLoc')));
if ~exist(out_dir,'dir') mkdir(out_dir); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% prep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% organize data

% load V first to avoid overwriting sortedResp
load(sprintf('uniLoc_sub%i_ses%s', sub, '-V'))
load(sprintf('uniLoc_sub%i_ses%s', sub, '-A'))

% set up
nRep = ExpInfo.nRep;
targInds = unique([sortedResp.target_idx]);
targDegs = unique([sortedResp.target_cm] * ScreenInfo.numPixels_perCM);
targNum = length(targInds);
rA = NaN(targNum,nRep);
rV1 = NaN(targNum,nRep);
rV2 = NaN(targNum,nRep);

for i = 1:targNum
    targInd = targInds(i);
    estDegsA = [sortedResp.response_pixel];
    estDegsV1 = [sortedReli1Resp.response_pixel];
    estDegsV2 = [sortedReli2Resp.response_pixel];

    rA(i,:) = estDegsA([sortedResp.target_idx] == targInd);
    rV1(i,:) = estDegsV1([sortedReli1Resp.target_idx] == targInd);
    rV2(i,:) = estDegsV2([sortedReli2Resp.target_idx] == targInd);
end

%% fit a linear line to both
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

coeffLabels.subjects = unique([coeffLabels.subjects, sub]);
coeffLabels.conditions = {'A', 'V1', 'V2'};
coeffLabels.Cofficients = {'Slope', 'Intercept'};
Bias.coeffLabels = coeffLabels;

Bias.coeff(sub, 1, :) = coefsA;
Bias.coeff(sub, 2, :) = coefsV1;
Bias.coeff(sub, 3, :) = coefsV2;

fitLabels.subjects = unique([coeffLabels.subjects, sub]);
fitLabels.conditions = {'A', 'V1', 'V2'};
fitLabels.fit = {'Slope', 'Intercept'};
Bias.coeffLabels = fitLabels;

Bias.fit(sub, 1, :) = fitRA;
Bias.fit(sub, 2, :) = fitSV1;
Bias.fit(sub, 3, :) = fitSV2;

save(out_dir, "Bias")

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




