function d_btst = btst_allData(D_matching, D_pointing, D_unimodal,...
    D_bimodal_unityJdg, D_bimodal_locResp, nbtst)

global nBtst
nBtst = nbtst; d_btst = cell(1,nbtst);
%------------------------------data structure------------------------------
%D_matching: cell, 1 x 4(psychometric curves)
%D_pointing: mat,  2 (1st: stimulus location; 2nd: resp) x 240 (trial number)
%D_unimodal: mat,  2 (1st: stimulus location; 2nd: resp) x 120 (trial number) 
%                  x 2 (modality, 1st: A, 2nd: V)
%D_bimodal_unityJdg, 2 (cond) x 2 (phase) x 4 (A locs) x 4 (V locs) x 20 (trial number)
%D_bimodal_locResp,  2 (cond) x 2 (phase) x 4 (A locs) x 4 (V locs) 
%                    x 2 (modality, 1st: A, 2nd: V) x 10 (trial number)
%--------------------------------------------------------------------------

%bootstrap data from the matching task
d_matching_btst = btst_matching(D_matching); 

%bootstrap data from the unimodal task
[d_pointing_btst,sigma_r_btst] = btst_pointing(D_pointing); 

%bootstrap data from the unimodal spatial-localization task
[d_unimodalA_btst, d_unimodalV_btst] = btst_unimodal(D_unimodal); 

%bootstrap data from the bimodal spatial-discrimination task
[d_unityJdg_btst,d_locResp_btst] = btst_unityJdg_locResp(D_bimodal_unityJdg,...
    D_bimodal_locResp); 

%put all bootstrapped data together in a cell
for n = 1:nBtst
    %select a set of bootstrapped data
    D.matching        = d_matching_btst{n};
    D.pointing        = squeeze(d_pointing_btst(n,:,:));
    D.sigma_r         = sigma_r_btst(n);
    D.unimodal(:,:,1) = squeeze(d_unimodalA_btst(n,:,:));
    D.unimodal(:,:,2) = squeeze(d_unimodalV_btst(n,:,:));
    D.bimodal_unity   = squeeze(d_unityJdg_btst(n,:,:,:,:,:));
    D.bimodal_locResp = squeeze(d_locResp_btst(n,:,:,:,:,:,:));
    D.s_A                   = unique(squeeze(D_unimodal(1,:,1)));
    D.s_V                   = unique(squeeze(D_unimodal(1,:,2)));
    D.numUnityTrialsPerLoc  = 20;
    d_btst{n} = D; clear D;
end

%--------------------------------------------------------------------------
%                       Bootstrap functions
%--------------------------------------------------------------------------
function d_matching_btst = btst_matching(d_matching)
global nBtst
numPF           = length(d_matching);    %4 psychometric functions
nTT             = size(d_matching{1},2); %80 per psychometric functions
d_matching_btst = cell(1,nBtst);
for i = 1:nBtst
    %resample with replacement
    d_matching_btst{i} = cell(1,numPF);
    idx_selected = randi(nTT,[numPF, nTT]); %4 x 80
    for j = 1:numPF;d_matching_btst{i}{j} = d_matching{j}(:,idx_selected(j,:));end
end

function [d_pointing_btst, sigma_r] = btst_pointing(d_pointing)
global nBtst
s_unique        = unique(d_pointing(1,:));             %8 stimulus locations
nT_per          = size(d_pointing,2)/length(s_unique); %30 trials per loc
d_pointing_btst = NaN(nBtst, 2, size(d_pointing,2));   
%2nd dim: (1st: stimulus location; 2nd: selected localization responses)
for j = 1:length(s_unique)
    %select responses that correspond to stimulus location = s_unique(j)
    resp_selected = d_pointing(2,d_pointing(1,:)==s_unique(j));
    for i = 1:nBtst
        numRand = randi(nT_per, 1, nT_per); %30 random integers (<=30)
        d_pointing_btst(i,1, ((j-1)*nT_per+1):(j*nT_per)) = ones(1,nT_per).*s_unique(j);
        d_pointing_btst(i,2, ((j-1)*nT_per+1):(j*nT_per)) = resp_selected(numRand);
    end
end
%calculate sigma_r
sigma_r   = NaN(1,nBtst);
for i = 1:nBtst
    locError          = squeeze(d_pointing_btst(i,2,:) - d_pointing_btst(i,1,:));
    SD                = std(locError);
    bool_nonoutliers  = (locError./SD < 3 & locError./SD > -3);
    cursor_locResp_rm = locError(bool_nonoutliers);
    %We use MLE to calculate sigma_r, which is sqrt(SSE/N)
    %We can't just use function std.m because the denominator is N-1
    sigma_r(i)        = sqrt(sum((cursor_locResp_rm - mean(cursor_locResp_rm)).^2)/...
                            length(cursor_locResp_rm));
end

function [d_unimodalA_btst, d_unimodalV_btst] = btst_unimodal(d_unimodal)
global nBtst
%2 (1st: stimulus locations; 2nd: localization responses) x 120 (trial number)
d_unimodalA = squeeze(d_unimodal(:,:,1)); 
d_unimodalV = squeeze(d_unimodal(:,:,2));
A_locs      = unique(d_unimodalA(1,:)); %4 unique auditory locs
V_locs      = unique(d_unimodalV(1,:)); %4 unique visual locs
nTT_perM    = size(d_unimodalA,2);      %# total trials per modality = 120
nT_per      = nTT_perM/length(A_locs);  %#trials per modality per location = 30
[d_unimodalA_btst, d_unimodalV_btst] = deal(NaN(nBtst,2, nTT_perM));
%2nd dim: (1st: stimulus location; 2nd: selected localization responses)
for j = 1:length(A_locs)
    %select responses that correspond to stimulus location = A_locs(j)
    respA_selected = d_unimodalA(2, abs((d_unimodalA(1,:) - A_locs(j)))<1e-4);
    respV_selected = d_unimodalV(2, abs((d_unimodalV(1,:) - V_locs(j)))<1e-4);
    for i = 1:nBtst
        numRand = randi(nT_per, 1, nT_per); %30 random integers (<=30)
        d_unimodalA_btst(i,1, ((j-1)*nT_per+1):(j*nT_per)) = ones(1,nT_per).*A_locs(j);
        d_unimodalV_btst(i,1, ((j-1)*nT_per+1):(j*nT_per)) = ones(1,nT_per).*V_locs(j);
        d_unimodalA_btst(i,2, ((j-1)*nT_per+1):(j*nT_per)) = respA_selected(numRand);
        d_unimodalV_btst(i,2, ((j-1)*nT_per+1):(j*nT_per)) = respV_selected(numRand);
    end
end

function [d_unityJdg_btst,d_locResp_btst] = btst_unityJdg_locResp(d_unityJdg,...
    d_locResp)
global nBtst
%get dimensionalities
lenC   = size(d_locResp,1); lenP   = size(d_locResp,2);
A_loc  = size(d_locResp,3); V_loc  = size(d_locResp,4);
lenM   = size(d_locResp,5); nT_per = size(d_locResp,6);
%initialization
d_unityJdg_btst = NaN(nBtst, lenC, lenP, A_loc, V_loc, lenM*nT_per);
d_locResp_btst  = NaN(nBtst, lenC, lenP, A_loc, V_loc, lenM, nT_per);
%for each dimension
for i = 1:lenC
    for j = 1:lenP
        for k = 1:A_loc
            for l = 1:V_loc
                for m = 1:lenM
                    %first 10 trials: A; last 10 trials: V
                    idx_lb = (m-1)*nT_per+1; %=1 when m = 1; =11 when m = 2;
                    idx_ub = nT_per*m;       %=10 when m = 1; =20 when m = 2;
                    d_unityJdg_selected = squeeze(d_unityJdg(i,j,k,l,idx_lb:idx_ub))'; 
                    d_locResp_selected  = squeeze(d_locResp(i,j,k,l,m,:))';
                    for n = 1:nBtst
                        numRand = randi(nT_per, 1, nT_per); %10 random integers (<=10)
                        d_unityJdg_btst(n,i,j,k,l,idx_lb:idx_ub) =...
                                        d_unityJdg_selected(numRand);
                        d_locResp_btst(n,i,j,k,l,m,:)  = ...
                                        d_locResp_selected(numRand);
                    end
                end
            end
        end
    end
end

