%% trying the same thing on data
% num_s: stimulus location combination
% 3 decision strategies
% 2 modalities(1 = aud, 2 = vis)
% 2 visual reliabilities
% num_rep
clear; clc;

sub_slc = 4;
ses_slc = 1:3;
total_num_rep = 16;

[org_resp, org_conf, org_err,org_uni, ExpInfo] = org_data(sub_slc,ses_slc,total_num_rep);


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
for modality = 1:2

    curr_resp = NaN(total_num_rep, 4);
    figure
    t = tiledlayout(2, 1);
    title(t, [modal_str{modality} ' locolization'])
    xlabel(t, 'Confidence radius')
    ylabel(t, 'Count');
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    % Loop over each auditory inde
    for vis_reliability = 1:2 %1:numel(ds_conf)
        % plot set up
        nexttile
        % Loop over each auditory index
        for i = 1:4
            curr_resp(:,i) = squeeze(org_conf(i, i, modality, vis_reliability, :));

        end
        h = histogram(curr_resp, 100);
        h.FaceColor = repmat(0.5, 1, 3);
        h.EdgeColor = 'none';
        xlim([0,40])
    end
end  


%% analyze data

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

%% plot localization response
edges = linspace(-512,512,128);
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

                curr_resp = squeeze(org_resp(a, v, modality, rel, :));
                h = histogram(curr_resp, edges);
                h.FaceColor = clt(rel+1,:);
                h.EdgeColor = 'none';
                h.FaceAlpha = 0.5;


                
            end                
            xline(sA(a),'--','LineWidth',lw,'Color',clt(1,:))
                xline(sV(v),'--','LineWidth',lw,'Color',clt(2,:))

                xlim([-512, 512])                
            if (a == 1) && (v == 4)
                    legend('high rel est','low rel est','Aud Stim', 'Vis Stim','Location','northeast')
                end
            hold off
        end
    end
end

%% plot confidence as a function of discrepancy and reliability

for d = 1 %1:numel(ds_conf)

    figure; hold on

    t = tiledlayout(2, 1);
    title(t, [])
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

end

%%
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

            histogram(squeeze(diff_conf{1, diff}(modality, rel, :)),'BinWidth',5);
            xlim([0, 250])

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
   
end


%% metamer trials to congruent trials

figure;
set(gcf, 'Position',[10 10 1200 300]); hold on
t = tiledlayout(1,4);

all_con = [];
all_incon = [];

for i  = 1:4

    slc_cue = 1;
    slc_rel = 1;
    slc_resp = squeeze(org_resp(:,:,slc_cue,slc_rel,:));
    con_conf = org_conf(i,i,slc_cue, slc_rel,:);
    all_con = [all_con; con_conf];

    con_resp = slc_resp(i,i,:);
    i_mu = mean(con_resp);
    i_sd = std(con_resp,[],'all');

    lb = i_mu - 10;%0.5*i_sd;
    ub = i_mu + 10;%0.5*i_sd;

    withinRange = (slc_resp>= lb) & (slc_resp<= ub);
    linearIndices = find(withinRange);

    [I1, I2, I3] = ind2sub(size(slc_resp), linearIndices);

    incon_conf = [];

    for j = 1:length(linearIndices)

        if I1(j) == i && I2(j) == i
        else

            incon_conf = [incon_conf; org_conf(I1(j), I2(j), slc_cue, slc_rel, I3(j))];
            all_incon = [all_incon; incon_conf];

        end

    end

    nexttile
    hold on
    h= histogram(con_conf,'BinWidth',5);
    h.FaceColor = clt(2,:);
    h.EdgeColor = 'none';
    h.FaceAlpha = 0.5;

    h= histogram(incon_conf,'BinWidth',5);
    h.FaceColor = repmat(0.5, 1, 3);
    h.EdgeColor = 'none';
    h.FaceAlpha = 0.5;

    xlim([0, 300])
    title(sprintf('estimate: %.2f - %.2f', lb, ub))

    if i == 4
        legend('congruent', 'incongruent')
    end
end


figure;
hold on
h= histogram(all_con,'BinWidth',5);
h.FaceColor = clt(2,:);
h.EdgeColor = 'none';
h.FaceAlpha = 0.5;

h= histogram(all_incon,'BinWidth',5);
h.FaceColor = repmat(0.5, 1, 3);
h.EdgeColor = 'none';
h.FaceAlpha = 0.5;
xlim([0, 300])
