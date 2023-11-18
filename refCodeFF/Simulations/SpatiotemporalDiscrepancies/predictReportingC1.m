clear all; close all; clc
%first define parameters
param.sigma_deltaT       = 80;
param.sigma_spatial_AV_A = 5;
param.sigma_spatial_AV_V = 1;
param.sigmaP_spatial     = 100;
param.muP_spatial        = 0;
param.sigmaP_deltaT_C1   = 100;
param.sigmaP_deltaT_C2   = 1e4;
param.muP_temporal       = 0;
pCommon      = 0.5;
spatialD     = -24:1:24;
temporalD    = -500:10:500;

combD        = combvec(spatialD, temporalD);
numTrials  = 100;
ds_unityJdg  = 'Euclidean distance';
c_spatiotemporal = 15; 
f_s = 1./5;
f_t = 1./500;
propC1 = NaN(length(spatialD), length(temporalD));

for i = 1:length(spatialD)
    for j = 1:length(temporalD)
        propC1(i,j) = predict_unityJdg_heuristic(param, numTrials, spatialD(i), temporalD(j),...
            pCommon, ds_unityJdg, NaN, NaN, c_spatiotemporal, f_s, f_t);
    end
end

%% plot
figure;
imagesc(temporalD, spatialD, propC1); colorbar
% yticks(linspace(1,length(spatialD),5)); 
% xticks(linspace(1,length(temporalD),5));
xlabel('Temporal offset'); ylabel('Spatial offset');



