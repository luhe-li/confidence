clear; close all;

%% out dir
cur_dir                          = pwd;
[project_dir, ~]                 = fileparts(cur_dir);
out_dir                          = fullfile(cur_dir, 's1Fig');
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% experimental info

% fix parameters for the experiment info
speaker_span       = 65.5 * 2; % cm
sitting_dist       = 113; % cm
screen_width       = 170; % cm
screenX            = 1024; % pixel
screen_mid         = screenX ./2;
num_speaker_int    = 15; % 15 intervals between 16 speakers
cm_per_aud_ind     = speaker_span / num_speaker_int;
pixel_per_cm       = screenX / screen_width;
aud_level          = [5 7 10 12];
fixP.screenX       = 1024;
fixP.x             = 1:1024; % the screen pixel space
center_x           = fixP.x - 512;

sA                 = round((aud_level - 8.5) .* cm_per_aud_ind .* pixel_per_cm + screen_mid); % in pixel. 8.5 being middle point
sV                 = sA;
aA                 = 1;
sAV                = combvec(sA, sV);
num_s              = size(sAV,2); % number of conditions

% fix parameters for the cost function
fixP.maxPoint     = 1;
fixP.minPoint     = 0.01;
fixP.elbow        = 256; % 1/4 of ScreenX is what we are using right now

% experimental info for analysis and plot
aud_locs          = sA - 512;
vis_locs          = sV - 512;
diffs       = zeros(length(aud_locs), length(vis_locs));
for i = 1:length(aud_locs)
    for j = 1:length(vis_locs)
        diffs(i, j) = aud_locs(i) - aud_locs(j);
    end
end
disc_locs   = unique(abs(diffs));

%% free parameters

bA                 = 0;
sigA               = 100;
sigVs              = [50,125];
muP                = screen_mid;
sigP               = 10000; % arbitarily using screen width here
pCommon            = 0.57; % only 1/4 of the trials are common cause so I assume this here

% derivative info
num_rep            = 1000;
ds_loc             = {'Model averaging E[V]','Model averaging MAP', 'Model selection', ...
    'Probability matching'};
ds_conf            = {'Model averaging', 'Model selection', ...
    'Probability matching_{select}', 'Probability matching_{mix}'};
cue_label   = {'Post-cue: A','Post-cue: V'};
num_cue            = numel(cue_label);
rel_label   = {'High reliability','Low reliability'};

%% simulate data

loc                = NaN(num_s,numel(ds_loc), num_cue, numel(sigVs), num_rep);
% num_s: stimulus location combination
% 3 decision strategies
% 2 modalities(1 = aud, 2 = vis)
% 2 visual reliabilities
% num_rep
conf               = NaN(num_s,numel(ds_conf), num_cue, numel(sigVs), num_rep);
unity              = NaN(num_s, numel(sigVs), num_rep);

for i = 1:num_s
    for v = 1:numel(sigVs)
        [loc(i,:,:,v,:), conf(i,:,:,v,:), unity(i,v,:), pdf(i,v)] = sim_loc_conf_unity_resp(pCommon, num_rep, sAV(1,i),...
            sAV(2,i), aA, bA, sigA, sigVs(v), muP, sigP, fixP);
    end
end

loc = loc - 512;


%% Plot set up

lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    128, 128, 128]./255; % light red

%% analyze data

allDiffs = unique(abs(sA' - sV));
num_diff = numel(allDiffs);

org_resp = NaN(numel(ds_loc),numel(sA), numel(sV), num_cue, numel(sigVs), num_rep);
org_conf = NaN(numel(ds_conf),numel(sA), numel(sV), num_cue, numel(sigVs), num_rep);
[mean_conf, sd_conf] = deal(NaN(numel(ds_conf),num_diff, num_cue, numel(sigVs)));

for d = 1:numel(ds_conf)

    % organize data

    if d == numel(ds_conf)
        d_loc = 3;
    else
        d_loc = d;
    end

    d_resp = squeeze(loc(:,d_loc,:,:,:));
    org_resp(d_loc,:,:,:,:,:) = reshape(d_resp, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

    d_conf = squeeze(conf(:,d,:,:,:));
    org_conf(d,:,:,:,:,:) = reshape(d_conf, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

    % organize response by discrepancy

    for i = 1:length(allDiffs)

        diff = allDiffs(i);

        % Find pairs of audIdx and visIdx that match this difference
        [audPairs, visPairs] = find(abs(sA' - sV) == diff);
        tempData = [];

        % For each pair, extract and store the corresponding data
        for j = 1:numel(audPairs)
            % Extract data for this specific audIdx and visIdx pair across all other dimensions
            tempData = cat(3, tempData, squeeze(org_conf(d, audPairs(j), visPairs(j), :, :, :)));
        end

        % store conf by discrepancy
        diff_conf{d,i} = tempData;

        % take average and sd
        mean_conf(d,i,:,:) = squeeze(mean(tempData,3));
        sd_conf(d,i,:,:) = squeeze(std(tempData,[],3));

    end

end

%% plot localization response

for d = 1:numel(ds_loc)

    % plot set up
    figure
    set(gcf, 'Position', get(0, 'Screensize')); hold on
    t = tiledlayout(4, 4);
    title(t, ['Auditory response, ' ds_loc{d}])
    xlabel(t, 'Stimulus Location (pixel)');
    ylabel(t, 'Count');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    % Loop over each auditory index
    for a = 1:length(sA)
        % Loop over each visual index
        for v = 1:length(sV)

            nexttile
            hold on

            for rel = 1:2

                curr_resp = squeeze(org_resp(d, a, v, 1, 1, :));
                h = histogram(curr_resp, 100);
                h.FaceColor = clt(rel+1,:);
                h.EdgeColor = 'none';
                h.FaceAlpha = 0.5;

                xline(aud_locs(a),'LineWidth',lw,'Color',clt(1,:))
                xline(vis_locs(v),'LineWidth',lw,'Color',clt(2,:))

                xlim([-512, 512])
                if (a == 1) && (v == 4)
                    legend('s_{hat}','Aud Stim', 'Vis Stim','Location','northeast')
                end

            end
        end
    end
end

%% plot confidence as a function of discrepancy and reliability

for d = 1:numel(ds_conf)

    figure; hold on

    t = tiledlayout(2, 1);
    title(t, ds_conf{d})
    xlabel(t, 'Audiovisual discrepancy (pixel)');
    ylabel(t, 'Confidence (pixel)');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    for cue = 1:num_cue
        nexttile
        title(cue_label{cue})
        hold on
        for rel = 1:numel(sigVs)
            errorbar(disc_locs, squeeze(mean_conf(d,:,cue,rel)), ...
                squeeze(sd_conf(d,:,cue,rel)),'Color',clt(rel+1,:),'LineWidth',1)
        end
        xticks(disc_locs)
    end
    legend(rel_label, 'Location', 'eastoutside');

    saveas(gca, fullfile(out_dir, sprintf('conf_%s',  ds_conf{d})), 'png')

end

%% plot confidence histogram as a function of discrepancy, fixing reliability

for d = 1:numel(ds_conf)

    figure;
    set(gcf, 'Position',[10 10 1200 400]); hold on
    t = tiledlayout(2, 5);
    title(t, ds_conf{d})
%     xlabel(t, 'Audiovisual discrepancy (pixel)');
    ylabel(t, 'Count');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    for rel = 1:2

        for diff = 1:num_diff
        
            nexttile
            hold on

            histogram(squeeze(diff_conf{d, diff}(cue, rel, :)),'BinWidth',5);
            xlim([0, 200])

            if diff == 1 && rel == 1
                ylabel('High visual reliability')
            elseif diff == 1 && rel == 2
                ylabel('Low visual reliability')
            end

            if rel == 2
                xlabel(sprintf('Discrepancy = %i', disc_locs(diff)))
            end

        end
    end

    saveas(gca, fullfile(out_dir, sprintf('hist_conf_aud_%s',  ds_conf{d})), 'png')
   
end


%% plot posterior pdf by conditions and strategies

org_pdf = reshape(pdf, [4, 4, 2]);


for d = 1:2%:numel(ds_conf)

    % plot set up
    figure
    set(gcf, 'Position', get(0, 'Screensize')); hold on
    t = tiledlayout(4, 4);
    title(t, ['Posterior of auditory stimulus location: ' ds_conf{d}])
    xlabel(t, 'Stimulus Location (pixel)');
    ylabel(t, 'Probability');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    % Loop over each auditory index
    for a = 1:length(sA)
        % Loop over each visual index
        for v = 1:length(sV)

            nexttile
            hold on

            for r = 1%:2

                if d == 1
                    post = org_pdf(a,v,r).MA_A;
                elseif d == 2
                    post = org_pdf(a,v,r).MS_A;
                elseif d == 3
                    post = org_pdf(a,v,r).PM_select_A;
                elseif d == 4
                    post = org_pdf(a,v,r).PM_match_A;
                end

                post = post(1,:);
                plot(center_x, post,'LineWidth', 1, 'Color', clt(1,:));

                xline(aud_locs(a),'LineWidth',lw,'Color',clt(1,:))
                xline(vis_locs(v),'LineWidth',lw,'Color',clt(2,:))
                xline(squeeze(org_resp(1, a, v, 1, r, 1)),'LineWidth',lw,'Color',clt(3,:))
                
                xlim([min(center_x), max(center_x)])

            end



            if (a == 1) && (v == 4)
                 legend('Posterior','Aud Stim', 'Vis Stim','estX','Location','northeast')
            end

        end
    end

    saveas(gca, fullfile(out_dir, sprintf('pdf_aud_%s',  ds_conf{d})), 'png')

end



