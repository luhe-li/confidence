%% trying the same thing on data
% num_s: stimulus location combination
% 3 decision strategies
% 2 modalities(1 = aud, 2 = vis)
% 2 visual reliabilities
% num_rep
clear; clc; 
sub_slc     = 4;
ses_slc     = 1;

% manage path
cur_dir      = pwd;
[project_dir, ~]= fileparts(fileparts(cur_dir));
out_dir      = fullfile(cur_dir, 's1Fig');
data_dir     = fullfile(project_dir, 'confidence','data','biLoc');
if ~exist(out_dir,'dir') mkdir(out_dir); end

% organize data
flnm        = sprintf('biLoc_sub%i_ses%i', sub_slc, ses_slc);
load(fullfile(data_dir, flnm))

% experiment info
aud_locs    = ExpInfo.speakerLocPixel(ExpInfo.audIdx);
vis_locs    = round(ExpInfo.targetPixel);
diffs       = zeros(length(aud_locs), length(vis_locs));
for i = 1:length(aud_locs)
    for j = 1:length(vis_locs)
        diffs(i, j) = aud_locs(i) - aud_locs(j);
    end
end
disc_locs   = unique(abs(diffs));

% conditions
seq         = ExpInfo.randAVIdx;
audIdx      = ExpInfo.audIdx;
visIdx      = ExpInfo.visIdx;
cueIdx      = ExpInfo.cueIdx;
visReliIdx  = ExpInfo.visReliIdx;
num_rep     = ExpInfo.nRep;
disc        = abs(seq(1,:) - seq(2,:));
discIdx     = unique(disc);
cue_label   = {'Post-cue: A','Post-cue: V'};
rel_label   = {'High visual reliability','Low visual reliability'};

% data
target      = [Resp.target_cm] * ScreenInfo.numPixels_perCM; % convert target to pixel, center as 0
resp        = [Resp.response_pixel] - ScreenInfo.xmid; % rescale response with center as 0
err         = abs(resp - target);
conf        = [Resp.conf_radius_cm];

% organize response
[org_resp, org_conf, org_err] = deal(NaN(numel(audIdx), numel(visIdx), numel(cueIdx), numel(visReliIdx), num_rep));
for trial = 1:length(resp)
    % Find indices for each condition
    aIdx = find(audIdx == seq(1, trial));
    vIdx = find(visIdx == seq(2, trial));
    cIdx = find(cueIdx == seq(3, trial));
    rIdx = find(visReliIdx == seq(4, trial));
    
    % Find the repetition index for the current combination
    repMatrix = squeeze(org_resp(aIdx, vIdx, cIdx, rIdx, :)); % Squeeze to remove singleton dimensions
    repIdx = find(isnan(repMatrix), 1); % Find the first NaN value
    
    % Insert the response
    org_resp(aIdx, vIdx, cIdx, rIdx, repIdx) = resp(trial);
    org_conf(aIdx, vIdx, cIdx, rIdx, repIdx) = conf(trial);
    org_err(aIdx, vIdx, cIdx, rIdx, repIdx) = err(trial);
end


% organize response by discrepancy
allDiffs = unique(abs(audIdx' - visIdx));

[mean_conf, sd_conf] = deal(NaN(numel(allDiffs), numel(cueIdx), numel(visReliIdx), num_rep));
for i = 1:length(allDiffs)
    diff = allDiffs(i);
    % Find pairs of audIdx and visIdx that match this difference
    [audPairs, visPairs] = find(abs(audIdx' - visIdx) == diff);
    tempData = NaN(numel(audPairs), numel(cueIdx), numel(visReliIdx), num_rep);
    % For each pair, extract and store the corresponding data
    for j = 1:numel(audPairs)
        % Extract data for this specific audIdx and visIdx pair across all other dimensions
        tempData(j,:,:,:) = squeeze(org_conf(audPairs(j), visPairs(j), :, :, :));
    end
    
    % Store the collected data in the main cell array
    mean_conf(i,:,:,:) = mean(tempData,1);
    sd_conf(i,:,:,:) = std(tempData,[],1);

end


%%
% Plot set up
nRep = 4;
lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    251, 154, 153]./255; % light red

curr_resp = NaN(nRep, length(aud_locs));
    figure
    t = tiledlayout(2, 1);
    title(t, 'Aud locolization' )
    xlabel(t, 'Confidence radius')
    ylabel(t, 'Count');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    % Loop over each auditory inde
for vis_reliability = 1:2 %1:numel(ds_conf)
    % plot set up
    nexttile

    title('Auditory response')
    xlabel('Stimulus Location (pixel)');
    ylabel('Count');

    % Loop over each auditory index
    for i = 1:length(aud_locs)
        curr_resp(:,i) = squeeze(org_conf(i, i, 1, vis_reliability, :));

    end
    h = histogram(curr_resp, 100);
    h.FaceColor = repmat(0.5, 1, 3);
    h.EdgeColor = 'none';
    xlim([0,40])

end



