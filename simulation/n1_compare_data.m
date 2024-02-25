clear; clc;
sub_slc = 4;
ses_slc = 1:3;
total_num_rep = 16;

[org_resp, org_conf, org_err,org_uni, ExpInfo] = org_data(sub_slc,ses_slc,total_num_rep);

org_uni_bool = org_uni>0;
%%
% Plot set up
nRep = 4;
lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    128 128 128]./255; % grey
modal_str = {"Audio","Visual"};
sA = ExpInfo.speakerLocPixel(ExpInfo.audIdx);
sV = sA;
allDiffs = unique(abs(sA' - sA));
sigVs = ExpInfo.visReliIdx;
cue_label   = {'Post-cue: A','Post-cue: V'};
num_cue            = numel(cue_label);
disc_locs   = unique(abs(allDiffs));
rel_label   = {'High reliability','Low reliability'};
num_diff = numel(allDiffs);


%%
[mean_conf, sd_conf] = deal(NaN(1,numel(allDiffs), num_cue, numel(sigVs)));

for d = 1 %1:numel(ds_conf)

    % organize response by discrepancy

    for i = 1:length(allDiffs)
        diff = allDiffs(i);
        % Find pairs of audIdx and visIdx that match this difference
        [audPairs, visPairs] = find(abs(sA' - sV) == diff);
        tempData = [];
        % For each pair, extract and store the corresponding data
        for j = 1:numel(audPairs)
            % Extract data for this specific audIdx and visIdx pair across all other dimensions
            tempData = cat(3, tempData, squeeze(org_conf(audPairs(j), visPairs(j), :, :, :)));
        end
        % store conf by discrepancy
        diff_conf{d,i} = tempData;
        % take average and sd
        mean_conf(d,i,:,:) = squeeze(mean(tempData,3));
        sd_conf(d,i,:,:) = squeeze(std(tempData,[],3));

    end

end

%%
mean_uni = NaN(1,numel(allDiffs), num_cue, numel(sigVs));

for d = 1 %1:numel(ds_conf)

    % organize response by discrepancy

    for i = 1:length(allDiffs)
        diff = allDiffs(i);
        % Find pairs of audIdx and visIdx that match this difference
        [audPairs, visPairs] = find(abs(sA' - sV) == diff);
        tempData = [];
        % For each pair, extract and store the corresponding data
        for j = 1:numel(audPairs)
            % Extract data for this specific audIdx and visIdx pair across all other dimensions
            tempData = cat(3, tempData, squeeze(org_uni_bool(audPairs(j), visPairs(j), :, :, :)));
        end
        % store conf by discrepancy
        diff_uni{d,i} = tempData;
        % take average and sd
        mean_uni(d,i,:,:) = squeeze(mean(tempData,3));


    end

end

%% plot uni boolean as a function of discrepancy and reliability

for d = 1 %1:numel(ds_conf)

    figure; hold on

    t = tiledlayout(2, 1);
    title(t, [])
    xlabel(t, 'Audiovisual discrepancy (pixel)');
    ylabel(t, '% Common Source');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    for cue = 1:num_cue
        nexttile
        title(cue_label{cue})
        hold on
        for rel = 1:numel(sigVs)
            plot(disc_locs, squeeze(mean_uni(d,:,cue,rel)), ...
                '-o','Color',clt(rel+1,:),'LineWidth',1)
        end
        xticks(disc_locs)
    end
    legend(rel_label, 'Location', 'eastoutside');

end

%%
% rel = 2;

for modality = 1:2 %:numel(ds_conf)

    figure;
    set(gcf, 'Position',[10 10 1200 400]); hold on
    t = tiledlayout(2, 5);
    title(t, [modal_str{modality} ' Confidence (pixel)'])
%     xlabel(t, 'Audiovisual discrepancy (pixel)');
    ylabel(t, 'Count');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    for rel = 1:2

        for diff = 1:num_diff
        
            nexttile
            hold on
            for scenario = 0:1
                curr_conf = squeeze(diff_conf{1, diff}(modality, rel, :));
                curr_uni_bool = squeeze(diff_uni{1, diff}(modality, rel, :));
                curr_conf(curr_uni_bool ~= scenario) = NaN;
                h = histogram(curr_conf,'BinWidth',5);
                h.FaceColor = clt(scenario+2,:);
                h.EdgeColor = 'none';
                h.FaceAlpha = 0.5;
            end
            xlim([0, 250])

            if diff == 1 && rel == 1
                ylabel('High Reli')
            elseif diff == 1 && rel == 2
                ylabel('Low Reli')
            end
            if diff == 5 && rel == 1
                legend('C=1','C=2')
            end
            if rel == 2
                xlabel(sprintf('Discrepancy = %i', disc_locs(diff)))
            end

        end
    end
   
end
%% plot localization response
edges = linspace(-512,512,128);
rel = 2;
for modality = 1:2
    % plot set up
    figure
    set(gcf, 'Position', get(0, 'Screensize')); hold on
    t = tiledlayout(4, 4);
    title(t, [modal_str{modality} ' locolization'])
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

            for scenario = 0:1

                curr_resp = squeeze(org_resp(a, v, modality, rel, :));
                curr_uni_bool = squeeze(org_uni_bool(a, v, modality, rel, :));
                curr_resp(curr_uni_bool ~= scenario) = NaN;
                h = histogram(curr_resp, edges);
                h.FaceColor = clt(scenario+2,:);
                h.EdgeColor = 'none';
                h.FaceAlpha = 0.5;


                
            end                
            xline(sA(a),'--','LineWidth',lw,'Color',clt(1,:))
                xline(sV(v),'--','LineWidth',lw,'Color',clt(2,:))

                xlim([-512, 512])                
            if (a == 1) && (v == 4)
                    legend('C = 1','C = 2','Aud Stim', 'Vis Stim','Location','northeast')
            end
            hold off
        end
    end
end

%% plot conf response
edges = linspace(0,200,40);
for modality = 1:2
    % plot set up
    figure
    set(gcf, 'Position', get(0, 'Screensize')); hold on
    t = tiledlayout(4, 4);
    title(t, [modal_str{modality} ' locolization'])
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

                curr_resp = squeeze(org_conf(a, v, modality, rel, :));
                h = histogram(curr_resp, edges);
                h.FaceColor = clt(rel+1,:);
                h.EdgeColor = 'none';
                h.FaceAlpha = 0.5;


                
            end                
            % xline(sA(a),'--','LineWidth',lw,'Color',clt(1,:))
            %     xline(sV(v),'--','LineWidth',lw,'Color',clt(2,:))
            % 
            %     xlim([-512, 512])                
            if (a == 1) && (v == 4)
                    legend('high rel est','low rel est','Aud Stim', 'Vis Stim','Location','northeast')
                end
            hold off
        end
    end
end


