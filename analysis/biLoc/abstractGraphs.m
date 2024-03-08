clear; clc; close all;
lw = 1;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [repmat(125, 1, 3); % blue
    227, 27, 27;  % dark red
    30, 120, 180]./255; % gray
%     251, 154, 153]./255; % light red
%%
cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(fileparts(cur_dir));

addpath(genpath(fullfile(project_dir,'simulation')))
addpath(genpath(fullfile(project_dir,'func')))

% out_dir                          = fullfile(cur_dir, 'n2Fig');
% if ~exist(out_dir,'dir') mkdir(out_dir); end


% fix parameters for the experiment info
speaker_span       = 65.5 * 2; % cm
sitting_dist       = 113; % cm
screen_width       = 170; % cm
screenX            = 1024; % pixel
screenXdeg         = rad2deg(atan(screen_width / 2 / sitting_dist)) .* 2;
screen_mid         = screenX ./2;
num_speaker_int    = 15; % 15 intervals between 16 speakers
cm_per_aud_ind     = speaker_span / num_speaker_int;
pixel_per_cm       = screenX / screen_width;
aud_level          = [6 8 9 11];
aud_VA             = -30:4:30;
deg_per_px         = screenXdeg / screenX;

fixP.screenX       = screenXdeg;
fixP.x             = -screenXdeg /2 : deg_per_px : screenXdeg/2; % the screen deg space
center_x           = 0;

sA                 = aud_VA(aud_level); % in deg. 8.5 being middle point
sV                 = sA;
aA                 = 1;
sAV                = combvec(sA, sV);
num_s              = size(sAV,2); % number of conditions

% fix parameters for the cost function
% fixP.maxPoint     = 1;
% fixP.minPoint     = 0.01;
% fixP.elbow        = 1000; % 1/4 of ScreenX is what we are using right now

% experimental info for analysis and plot
aud_locs          = sA;
vis_locs          = sV;
diffs             = zeros(length(aud_locs), length(vis_locs));
for i = 1:length(aud_locs)
    for j = 1:length(vis_locs)
        diffs(i, j) = aud_locs(i) - vis_locs(j);
    end
end
disc_locs   = unique(abs(diffs));


%% simulation
figure;
set(gcf, 'Position', [0 0 800 400]);
hold on

t = tiledlayout(2, 3);

xlabel(t, 'Audiovisual spatial discrepancy (deg)');
ylabel(t, 'Proportion of reporting confident');
title(t,'Model Simulation')
t.TileSpacing = 'compact';
t.Padding = 'compact';
tileIndices = [1 2 3; 4 5 6];
model_label = {"Heuristic","Suboptimal","Optimal"};
for model_ind = 1:numel(model_label)

bA                 = 0;
sigA               = [1.2;4;1];
sigVs              = [0.7,1;1,10;0.5,0.7];
muP                = 0;
sigP               = 10000; % arbitarily using screen width here
pCommon            = [0.57;0.01;0.57]; % only 1/4 of the trials are common cause so I assume this here
criterion = [1.45,0.5;13.8,13;1.05,0.3]; % eyeballed arbitary criterion 
lapse = 0.05;
indVar = "p"; 
% derivative info
num_rep            = 100;
ds_loc             = {'Model averaging E[V]','Model averaging MAP', 'Model selection','Probability Matching'};
ds_conf            = {'M1','M2','M3'};
cue_label   = {'Auditory confidence','Visual confidence'};
rel_label   = {'High vis reliability','Low vis reliability'};
num_cue            = numel(cue_label);

loc                = NaN(num_s,numel(ds_loc), num_cue, numel(rel_label), num_rep);
% num_s: stimulus location combination
% 3 decision strategies
% 2 modalities(1 = aud, 2 = vis)
% 2 visual reliabilities
% num_rep

variance           = NaN(num_s,numel(ds_conf), num_cue, numel(rel_label), num_rep);
unity              = NaN(num_s, numel(rel_label), num_rep);
p_conf             = NaN(num_s,numel(ds_conf), num_cue, numel(rel_label)); % this is summed per condition so just one less dimension
for i = 1:num_s
    for v = 1:numel(rel_label)
        [loc(i,:,:,v,:), unity(i,v,:), pdf(i,v), variance(i,:,:,v,:), p_conf(i,:,:,v)] = sim_loc_conf_unity_resp_EVResp(pCommon(model_ind), num_rep, sAV(1,i),...
            sAV(2,i), aA, bA, sigA(model_ind), sigVs(model_ind,v), muP, sigP, fixP, criterion, lapse);
    end
end

allDiffs = unique(abs(sA' - sV));
num_diff = numel(allDiffs);

for d = 1:numel(ds_conf)

    % organize data
    d_p_conf = squeeze(p_conf(:,d,:,:));
    org_p_conf(d,:,:,:,:) = reshape(d_p_conf, [numel(sA), numel(sV), num_cue, numel(rel_label)]);

    % organize response by discrepancy

    for i = 1:length(allDiffs)
        diff = allDiffs(i);
        % Find pairs of audIdx and visIdx that match this difference
        [audPairs, visPairs] = find(abs(sA' - sV) == diff);
        tempData_p_conf = [];
        % For each pair, extract and store the corresponding data
        for j = 1:numel(audPairs)
            tempData_p_conf = cat(3, tempData_p_conf, squeeze(org_p_conf(d, audPairs(j), visPairs(j), :, :)));
        end

        diff_p_conf{d,i} = tempData_p_conf;
        mean_p_conf(d,i,:,:) = squeeze(mean(tempData_p_conf,3));
        sd_p_conf(d,i,:,:) = squeeze(std(tempData_p_conf,[],3));
    end

end

    mean_iv = mean_p_conf;
    sd_iv = sd_p_conf;
    indVarName = "P(report Confident)";

    for cue = 1:num_cue
        nexttile(tileIndices(cue,model_ind))
        hold on
        if model_ind == 1
            ylabel(cue_label{cue})
        end
        if cue == 1
            title([model_label{model_ind} ' model'])
        end
        for rel = 1:numel(rel_label)
                plot(disc_locs, squeeze(mean_iv(model_ind,:,cue,rel)),'-o', ...
                    'Color',[clt(rel+1,:), 0.4],'LineWidth',lw,'DisplayName',rel_label{rel}, ...
                    'MarkerFaceColor',clt(rel+1,:),'MarkerSize',3)
                ylim([0,1])
            
        end
        xticks(round(disc_locs))
        
    end
    % legend(rel_label, 'Location', 'eastoutside');

    % saveas(gca, fullfile(out_dir, sprintf('conf_%s',  ds_conf{d})), 'png')

end
lgd = legend();
lgd.Layout.Tile = 'east';
%% data
ses_slc_struct(1).ses_slc = 1:3;
ses_slc_struct(2).ses_slc = 1:2;
ses_slc_struct(3).ses_slc = 1:2;
sub_slc_vector = [2,6,4];

figure; 
set(gcf, 'Position', [0 0 800 400]);
hold on 

t = tiledlayout(2, 3);

xlabel(t, 'Audiovisual spatial discrepancy (deg)');
ylabel(t, 'Proportion of reporting confident');
title(t,'Preliminary Data')
t.TileSpacing = 'compact';
t.Padding = 'compact';
tileIndices = [1 2 3; 4 5 6];
model_label = {"Heuristic","Suboptimal","Optimal"};
for model_ind = 1:3
    sub_slc = sub_slc_vector(model_ind);
    ses_slc     = [ses_slc_struct(model_ind).ses_slc];

% manage path
cur_dir      = pwd;
[project_dir, ~]= fileparts(fileparts(cur_dir));
% out_dir      = fullfile(cur_dir, 's1Fig');
addpath(genpath(fullfile(project_dir,'func')))
data_dir     = fullfile(project_dir, 'data','biLoc');
% if ~exist(out_dir,'dir') mkdir(out_dir); end

% organize data
[bi_resp, bi_conf, bi_err, ExpInfo] = org_data(sub_slc,ses_slc,'biLoc');
[uni_resp, uni_conf, uni_err, uniExpInfo, ~, ScreenInfo] = org_data(sub_slc,[],'uniLoc');

% conditions
seq         = ExpInfo.randAVIdx;
audIdx      = ExpInfo.audIdx;
visIdx      = ExpInfo.visIdx;
cueIdx      = ExpInfo.cueIdx;
visReliIdx  = ExpInfo.visReliIdx;
num_rep     = ExpInfo.nRep;
deg_per_px  = rad2deg(atan(170 / 2 / ExpInfo.sittingDistance)) .* 2 / ScreenInfo.xaxis;
% cue_label   = {'Auditory Confidence','Visual Confidence'};
% rel_label   = {'High vis','Low vis'};
aud_locs    = ExpInfo.speakerLocVA(ExpInfo.audIdx);
remapped_vis_locs = ExpInfo.targetPixel .* deg_per_px;

% summarize unimodal data
% condition (A,V1,V2) x loc (4) x rep
uni_pconf = mean(mean(uni_conf,3),2);
uni_idx = [1,1;2,3];

% organize by discrepancy
[conf_by_diff, diff_locs] = org_by_diffs(bi_conf, aud_locs); % {diff} cue x reliability x rep 

% summarize ventriloquist effect
% condition (A,V1,V2) x loc (4) x rep
a_uni = mean(squeeze(uni_resp(1,:,:)),2);
v1_uni = mean(squeeze(uni_resp(2,:,:)),2);
v2_uni = mean(squeeze(uni_resp(3,:,:)),2);
uni = {a_uni, a_uni; v1_uni', v2_uni'};
for cue = 1:2
    for rel = 1:2
        subset_bi_resp = squeeze(mean(bi_resp(:, :, cue, rel, :),5));
        if cue == 1
            uni_repeated = repmat(uni{cue, rel}, [1, size(bi_resp, 2)]);
            % ve(:,:,cue, rel) = subset_bi_resp - aud_locs';
        elseif cue == 2
            % uni_repeated = repmat(uni{cue, rel}, [size(bi_resp, 1),1]);
            uni_repeated = ExpInfo.targetPixel .* deg_per_px;
            % ve(:,:,cue, rel) = subset_bi_resp - aud_locs;
        end
        ve(:,:,cue, rel) = (subset_bi_resp - uni_repeated);
    end
end
[ve_by_diff, raw_diff] = org_by_raw_diffs(ve, aud_locs); % {diff} cue x reliability x rep
raw_diff(2,:) = -raw_diff;
raw_diff = -raw_diff;
for cue = 1:numel(cueIdx)
    nexttile(tileIndices(cue,model_ind))
    % title(cue_label{cue})
    if model_ind == 1
        ylabel(cue_label{cue})
    end  
    if cue == 1
        title([model_label{model_ind} ' participant'])
    end
    hold on
    for rel = 1:numel(visReliIdx)
        [mean_p, sd_p] = deal(NaN(1, numel(diff_locs)));
        for diff = 1:numel(diff_locs)
            est = squeeze(conf_by_diff{diff}(cue, rel, :))';
            p = sum(est)/numel(est);
            mean_p(diff) = p;
            sd_p(diff) = sqrt((p*(1-p))/numel(est));
        end
        plot(diff_locs, mean_p,'-o', 'Color',[clt(rel+1,:), 0.4],'LineWidth',lw,'DisplayName',rel_label{rel},'MarkerFaceColor',clt(rel+1,:),'MarkerSize',3)
        patch([diff_locs, fliplr(diff_locs)], ...
            [mean_p - sd_p, fliplr(mean_p + sd_p)], ...
            clt(rel+1,:),'EdgeColor','none','FaceAlpha',0.1, ...
            'HandleVisibility', 'off')
        
        if cue == 1
            yline(uni_pconf(1),'--','Color',clt(1,:),'LineWidth',lw,'DisplayName','Aud unimodal')
            
        else
            yline(uni_pconf(rel+1),'--','Color',clt(rel+1,:),'LineWidth',lw,'DisplayName',[rel_label{rel}])
        end
        
        ylim([0, 1])
        xlim([min(diff_locs), max(diff_locs)])
    end
    xticks(round(diff_locs))
    % legend({'High visual reliability','Unimodal baseline','Low visual reliability'}, 'Location', 'eastoutside');
end

end

lgd = legend();
lgd.Layout.Tile = 'east';
