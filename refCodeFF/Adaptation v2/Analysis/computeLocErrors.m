function [locErrors_A, locErrors_A_SD, locErrors_A_SD_btst,locErrors_A_SD_btst_CI,...
    locErrors_V, locErrors_V_SD, locErrors_V_SD_btst, locErrors_V_SD_btst_CI] =...
    computeLocErrors(sA, sV, r, m, rA_uni, rV_uni, nBtst, trials_slc)

if nargin < 7 || isempty(nBtst); nBtst = 1e3; end 
if nargin < 8; trials_slc = 'all'; end

%auditory/visual stimulus locations
sA_slc    = sA(m == 1); %80 trials
sV_slc    = sV(m == 2); 
%auditory/visual localization responses
rA        = r(m == 1);
rV        = r(m == 2);
nT        = length(sA_slc);
if strcmp(trials_slc, '1stHalf')
    sA_slc = sA_slc(1:nT/2); sV_slc = sV_slc(1:nT/2);
    rA     = rA(1:nT/2);     rV     = rV(1:nT/2);
elseif strcmp(trials_slc, '2ndHalf')
    sA_slc = sA_slc((nT/2+1):end); sV_slc = sV_slc((nT/2+1):end);
    rA     = rA((nT/2+1):end);     rV     = rV((nT/2+1):end);
elseif strcmp(trials_slc, '1stQuarter')
    sA_slc = sA_slc(1:nT/4); sV_slc = sV_slc(1:nT/4);
    rA     = rA(1:nT/4);     rV     = rV(1:nT/4);
elseif strcmp(trials_slc, '2ndQuarter')
    sA_slc = sA_slc((nT/4+1):(nT/2)); sV_slc = sV_slc((nT/4+1):(nT/2));
    rA     = rA((nT/4+1):(nT/2));     rV     = rV((nT/4+1):(nT/2));    
elseif strcmp(trials_slc, '3rdQuarter')
    sA_slc = sA_slc((nT/2+1):(3*nT/4)); sV_slc = sV_slc((nT/2+1):(3*nT/4));
    rA     = rA((nT/2+1):(3*nT/4));     rV     = rV((nT/2+1):(3*nT/4));
elseif strcmp(trials_slc, '4thQuarter')
    sA_slc = sA_slc((3*nT/4+1):end); sV_slc = sV_slc((3*nT/4+1):end);
    rA     = rA((3*nT/4+1):end);     rV     = rV((3*nT/4+1):end);    
end

%In congruent condition, there are 4 unique visual and 4 unique
%auditory locations; in incongruent condition, there are 8 unique
%visual locations (because they were chosen to be around the auditory
%stimuli)
find_nT   = @(v, s) sum(s == v);
sA_unique = unique(sA_slc); 
nT_perA   = arrayfun(@(ss) find_nT(sA_slc, sA_unique(ss)), 1:length(sA_unique));
idx_A     = [0, cumsum(nT_perA)];
sV_unique = unique(sV_slc); 
nT_perV   = arrayfun(@(ss) find_nT(sV_slc, sV_unique(ss)), 1:length(sV_unique));
idx_V     = [0, cumsum(nT_perV)];
%compute the localization errors (bimodal auditory localization
%responses - mean unimodal localization responses)
locErrors_A_all = NaN(1,length(sA_slc));
locErrors_A_all_btst = NaN(nBtst, length(sA_slc));
for i = 1:length(sA_unique)
    deMeanResp_A = rA(sA_slc == sA_unique(i)) - rA_uni(i);
    locErrors_A_all((idx_A(i)+1):idx_A(i+1))  = deMeanResp_A;
    %bootstrap data
    for j = 1:nBtst
        btst_idx = randi(nT_perA(i), [1, nT_perA(i)]);
        locErrors_A_all_btst(j, (idx_A(i)+1):idx_A(i+1)) = deMeanResp_A(btst_idx);
    end
end

locErrors_A         = mean(locErrors_A_all); 
locErrors_A_SD      = sqrt(sum((locErrors_A_all - locErrors_A).^2)/...
                        length(locErrors_A_all)); 
locErrors_A_btst    = mean(locErrors_A_all_btst,2);
locErrors_A_SD_btst = sqrt(sum((locErrors_A_all_btst - locErrors_A_btst).^2, 2)/...
                        size(locErrors_A_all_btst,2)); 
locErrors_A_SD_btst_CI = get95CI(locErrors_A_SD_btst, [0.025, 0.975], nBtst);

%compute the localization errors (bimodal visual localization
%responses - mean unimodal localization responses)
if length(sV_unique) == length(sA_unique); pred_rV_uni = rV_uni;
else
    estP = polyfit(rV_uni, -12:8:12, 1);
    pred_rV_uni = polyval(estP, sV_unique);
end

locErrors_V_all = NaN(1,length(sA_slc));
locErrors_V_all_btst = NaN(nBtst, length(sV_slc));
for i = 1:length(sV_unique)
    deMeanResp_V   = rV(sV_slc == sV_unique(i)) - pred_rV_uni(i);
    locErrors_V_all((idx_V(i)+1):idx_V(i+1)) = deMeanResp_V;

    %bootstrap data
    for j = 1:nBtst
        btst_idx = randi(nT_perV(i), [1, nT_perV(i)]);
        locErrors_V_all_btst(j,(idx_V(i)+1):idx_V(i+1)) = ...
            deMeanResp_V(btst_idx);
    end
end
locErrors_V         = mean(locErrors_V_all); 
locErrors_V_SD      = sqrt(sum((locErrors_V_all - mean(locErrors_V)).^2)/...
                        length(locErrors_V_all));    
locErrors_V_btst    = mean(locErrors_V_all_btst,2);
locErrors_V_SD_btst = sqrt(sum((locErrors_V_all_btst - locErrors_V_btst).^2, 2)/...
                        size(locErrors_V_all_btst,2)); 
locErrors_V_SD_btst_CI = get95CI(locErrors_V_SD_btst, [0.025, 0.975], nBtst);


function bds = get95CI(v, pct, btst)
v_sort = sort(v);
lb     = floor(btst*pct(1));
ub     = ceil(btst*pct(2));
bds    = [v_sort(lb), v_sort(ub)];

