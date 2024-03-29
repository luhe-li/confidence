clear; close all;


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
sigVs              = [25,150];
muP                = screen_mid;
sigP               = 10000; % arbitarily using screen width here
pCommon            = 0.57; % only 1/4 of the trials are common cause so I assume this here

% derivative info
num_rep            = 1000;
ds_loc             = {'Model averaging', 'Model selection', ...
    'Probability matching'};
ds_conf            = {'Model averaging', 'Model selection', ...
    'Probability matching_{select}', 'Probability matching_{mix}'};
cue_label   = {'Post-cue: A','Post-cue: V'};
num_cue            = numel(cue_label);

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
        [loc(i,:,:,v,:), conf(i,:,:,v,:), unity(i,v,:)] = sim_loc_conf_unity_resp(pCommon, num_rep, sAV(1,i),...
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
    251, 154, 153]./255; % light red

%% analyze data

for d = 1:numel(ds_loc)

    % organize data

    if d ~= numel(ds_conf)
        d_loc = d;
    else
    end

    d_resp = squeeze(loc(:,d_loc,:,:,:));
    org_resp = reshape(d_resp, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

    d_conf = squeeze(conf(:,d,:,:,:));
    org_conf = reshape(d_conf, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);


    % organize response by discrepancy
    allDiffs = unique(abs(sA' - sV));

    [mean_conf, sd_conf] = deal(NaN(numel(allDiffs), num_cue, numel(sigVs), num_rep));
    for i = 1:length(allDiffs)
        diff = allDiffs(i);
        % Find pairs of audIdx and visIdx that match this difference
        [audPairs, visPairs] = find(abs(sA' - sV) == diff);
        tempData = NaN(numel(audPairs), num_cue, numel(sigVs), num_rep);
        % For each pair, extract and store the corresponding data
        for j = 1:numel(audPairs)
            % Extract data for this specific audIdx and visIdx pair across all other dimensions
            tempData(j,:,:,:) = squeeze(org_conf(audPairs(j), visPairs(j), :, :, :));
        end

        % Store the collected data in the main cell array
        mean_conf(i,:,:,:) = mean(tempData,1);
        sd_conf(i,:,:,:) = std(tempData,[],1);

    end

end

%% plot localization response

for d = 1:numel(ds_loc)
    % plot set up
    figure
    set(gcf, 'Position', get(0, 'Screensize')); hold on
    t = tiledlayout(4, 4);
    title(t, ['Auditory response with ' ds_loc{d}])
    xlabel(t, 'Stimulus Location (pixel)');
    ylabel(t, 'Count');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    % Loop over each auditory index
    for a = 1:length(sA)
        % Loop over each visual index
        for v = 1:length(sV)
            curr_resp = squeeze(org_resp(a, v, 1, :, :));
            nexttile
            h = histogram(curr_resp, 100);
            h.FaceColor = repmat(0.5, 1, 3);
            h.EdgeColor = 'none';

            xline(aud_locs(a),'LineWidth',lw,'Color',clt(1,:))
            xline(vis_locs(v),'LineWidth',lw,'Color',clt(2,:))

            xlim([-512, 512])
            if (a == 1) && (v == 4)
                legend('s_{hat}','Aud Stim', 'Vis Stim','Location','northeast')
            end
        end
    end
end

%% plot confidence for congruent trials

for d = 1 %1:numel(ds_conf)

    % plot set up
    figure
    set(gcf, 'Position', get(0, 'Screensize')); hold on
    t = tiledlayout(4, 1);
    title(t, ['Auditory response with ' ds_conf{d}])
    xlabel(t, 'Stimulus Location (pixel)');
    ylabel(t, 'Count');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    % Loop over each auditory index
    for a = 1:length(sA)
        % Loop over each visual index
        for v = 1:length(sV)
            if a == v  % limiting only congurent trials
                curr_resp = squeeze(org_conf(a, v, 1, :, :));
                nexttile
                h = histogram(curr_resp, 100);
                h.FaceColor = repmat(0.5, 1, 3);
                h.EdgeColor = 'none';

                xline(aud_locs(a),'LineWidth',lw,'Color',clt(1,:))
                xline(vis_locs(v),'LineWidth',lw,'Color',clt(2,:))

                xlim([0, 512])
            end
        end
        % if (a == 1) && (v == 4)
        %     legend('s_{hat}','Aud Stim', 'Vis Stim','Location','northeast')
        % end
    end

end

%%
nRep = num_rep;
d = 1;
curr_resp = NaN(nRep, length(aud_locs));
    d_conf = squeeze(conf(:,d,:,:,:));
    org_conf = reshape(d_conf, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

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
    xlim([0,120])

end


%%

for d = 1 %1:numel(ds_loc)

    % organize data

    if d ~= numel(ds_conf)
        d_loc = d;
    else
    end

    d_resp = squeeze(loc(:,d_loc,:,:,:));
    org_resp = reshape(d_resp, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);

    d_conf = squeeze(conf(:,d,:,:,:));
    org_conf = reshape(d_conf, [numel(sA), numel(sV), num_cue, numel(sigVs), num_rep]);


    % organize response by discrepancy
    allDiffs = unique(abs(sA' - sV));

    [mean_conf, sd_conf] = deal(NaN(numel(allDiffs), num_cue, numel(sigVs), num_rep));
    for i = 1:length(allDiffs)
        diff = allDiffs(i);
        % Find pairs of audIdx and visIdx that match this difference
        [audPairs, visPairs] = find(abs(sA' - sV) == diff);
        tempData = NaN(numel(audPairs), num_cue, numel(sigVs), num_rep);
        % For each pair, extract and store the corresponding data
        for j = 1:numel(audPairs)
            % Extract data for this specific audIdx and visIdx pair across all other dimensions
            tempData(j,:,:,:) = squeeze(org_conf(audPairs(j), visPairs(j), :, :, :));
        end

        % Store the collected data in the main cell array
        mean_conf(i,:,:,:) = mean(tempData,1);
        sd_conf(i,:,:,:) = std(tempData,[],1);

    end

    for a = 1:length(aud_locs)
        % Loop over each visual index
        for v = 1:length(vis_locs)
            
                curr_resp(a,v,:,:) = squeeze(org_resp(a, v, 1, :, :));
                curr_conf(a,v,:,:) = squeeze(org_conf(a, v, 1, :, :));
                curr_err(a,v,:,:) = curr_resp(a,v,:,:) - aud_locs(a);
        end
        % if (a == 1) && (v == 4)
        %     legend('s_{hat}','Aud Stim', 'Vis Stim','Location','northeast')
        % end
    end
end



%%
figure
set(gca, 'LineWidth', lw, 'FontSize', fontSZ, 'TickDir', 'out')
set(gcf, 'Position',[0 0 500 400])
hold on 
scatter(err, conf,'filled','k')
[linearFit, gof] = fit(err', conf', 'poly1');
plot(linearFit, 'r-');
legend off
xlabel('Error (pixel)')
ylabel('Confidence (pixel)')









