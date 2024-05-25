clear; close all;

%% set recovery parameters

num_rep = 100; % experimental repetition per condition
num_run = 4; % run of fits

sampleGT = false; % use a fixed set of GT or sample num_sample from GT range
if sampleGT; num_sample = 100; else; num_sample = 1; end
checkFakeData = true;

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(fullfile(project_dir, 'bads')))
addpath(genpath(fullfile(project_dir, 'modelFit', 'biLocV3_4ratings')))

%% plot set up

lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;

%% experiment info

% fix parameters for the experiment info
speaker_span          = 65.5 * 2; % cm
sitting_dist          = 113; % cm
screen_width          = 170; % cm
screenX               = 1024; % pixel
screenXdeg            = rad2deg(atan(screen_width / 2 / sitting_dist)) .* 2;
screen_mid            = screenX ./2;
num_speaker_int       = 15; % 15 intervals between 16 speakers
aud_level             = [6 8 9 11];
aud_VA                = -30:4:30;
deg_per_px            = screenXdeg / screenX;

sA                    = aud_VA(aud_level); % in deg. 8.5 being middle point
sV                    = sA; % assume perceptually aligned
sAV                   = combvec(sA, sV);
num_s                 = size(sAV,2); % number of conditions

diff                 = zeros(length(sA), length(sV));
for ii                 = 1:length(sA)
    for j                 = 1:length(sV)
        diff(ii, j)           = sA(ii) - sV(j);
    end
end
abs_diff             = unique(abs(diff))';

%% model info

if ~sampleGT

    %        aA,    bA, sigV1, dsigA, dsigV2, sigP,   pCC,  sigC,    c1,  dc2,  dc3
    GT = {[   1,   0.1,     1,   1.2,    1.5,    8,  0.57,   0.3,   0.5, 1, 1],...% Heuristic
        [     1,   0.1,   0.5,   1.5,    1.8,    8,   0.8,   0.5,  0.48, 0.5, 0.5],...% Suboptimal
        [     1,   0.1,     1,   1.5,    1.8,    8,  0.57,   0.5,  0.48, 0.5, 0.5]}; % Optimal

else

    % TO-DO: sample from a GT range

end

ds_loc                = {'Model averaging'};
ds_conf               = {'Heuristic','Suboptimal','Optimal'};
num_model             = numel(ds_conf);
cue_label             = {'Post-cue: A','Post-cue: V'};
num_cue               = numel(cue_label);
rel_label             = {'High reliability','Low reliability'};

%% run

flnm = sprintf('recoveryResults_rep%i_sample%i_run%i', num_rep, num_sample, num_run);

if exist(flnm,'file')

    load(flnm);

else

    [saveData, saveModel, pred] = deal(cell(num_model, num_sample));

    for i = 1:num_sample

        for d = 3

            % specific GT for this sample and model
            i_gt = GT{i, d};

            % assign simulation parameters
            num_para              = length(i_gt);
            aA                    = i_gt(1);
            bA                    = i_gt(2);
            sigV1                 = i_gt(3);
            sigA                  = i_gt(4);
            sigV2                 = i_gt(5);
            sigVs                 = [sigV1, sigV2];
            sigP                  = i_gt(6);
            pCommon               = i_gt(7);
            sigC                  = i_gt(8);
            c1                    = i_gt(9);
            c2                    = i_gt(9) + i_gt(10);
            c3                    = i_gt(9) + i_gt(10) + i_gt(11);
            lapse                 = 0.02;
            muP                   = 0;

            % num_s: stimulus location combination
            % 3 confidence decision strategies
            % 2 modalities(1 = aud, 2 = vis)
            % 2 visual reliabilities
            % num_rep
            [loc, conf] = deal(NaN(num_s, num_cue, numel(sigVs), num_rep));

            for a = 1:numel(sA)

                for v = 1:numel(sV)

                    for r = 1:numel(sigVs)

                        fixP.sA = sA(a);
                        fixP.sV = sV(v);
                        fixP.model_ind = d;
                        fixP.sigMotor = 1.36; % in deg, emperical motor noise averaged from first four participants
                        fixP.num_rep = num_rep;

                        [org_loc(a,v,:,r,:), org_conf(a,v,:,r,:)] = simAllModels(...
                            aA, bA, sigA, sigVs(r), muP, sigP, pCommon, sigC, c1, c2, c3, lapse, fixP);

                    end

                end

            end

            %% check fake data

            if checkFakeData

                %% check localization data

                uni_loc = zeros(size(org_loc));

                loc_a = repmat(sA',[1,numel(sV)]);
                loc_v = repmat(sV,[numel(sA),1]);

                uni_loc(:,:,1,1:2,:) = repmat(loc_a, [1, 1, 1, 2, num_rep]);
                uni_loc(:,:,2,1:2,:) = repmat(loc_v, [1, 1, 1, 2, num_rep]);

                % loc at uni minus loc at bi
                ve =  mean(org_loc,5) - mean(uni_loc, 5);

                % diff x cue x reliability
                [ve_by_raw_diff, raw_diff] = org_by_raw_diffs_4D(ve, sA);

                % assume participants localized perfectly in the unisensory
                % condition
                figure; hold on
                t = tiledlayout(1, 2);
                title(t,sprintf('%s, rep: %i', ds_conf{d}, num_rep))
                xlabel(t, 'Audiovisual discrepancy (V-A, deg)');
                ylabel(t, 'Shift of localization');
                t.TileSpacing = 'compact';
                t.Padding = 'compact';

                for cue = 1:num_cue
                    nexttile
                    title(cue_label{cue})
                    axis equal
                    hold on

                    for rel = 1: numel(sigVs)

                        i_ve = squeeze(ve_by_raw_diff(:, cue, rel));
                        plot(raw_diff, i_ve, 'Color',clt(rel+1,:))

                    end
                    xlim([min(raw_diff), max(raw_diff)])
                    xticks(raw_diff)
                    yline(0,'--')
                    if cue == 1
                        plot(raw_diff, raw_diff,'k--')
                    else
                        plot(raw_diff, -raw_diff,'k--')
                    end
                end

                %% check confidence data

                % organize localization data: {diff} cue x reliability x rep
                [conf_by_diff, all_diffs] = org_by_diffs(org_conf, sA);

                figure; hold on
                t = tiledlayout(2, 1);
                title(t,sprintf('%s, rep: %i', ds_conf{d}, num_rep))
                xlabel(t, 'Absolute audiovisual discrepancy (deg)');
                ylabel(t, 'Proportion of confidence report');
                t.TileSpacing = 'compact';
                t.Padding = 'compact';

                for cue = 1:num_cue
                    nexttile
                    title(cue_label{cue})
                    hold on
                    for rel = 1: numel(sigVs)
                        [p_conf, se_conf] = deal(NaN(1, numel(abs_diff)));
                        for diff = 1:numel(abs_diff)
                            i_conf = squeeze(conf_by_diff{diff}(cue, rel, :))';
                            p = sum(i_conf)/(numel(i_conf));
                            p_conf(diff) = p;
                            se_conf(diff) = sqrt((p*(1-p))/numel(i_conf));
                        end
                        plot(abs_diff, p_conf, 'Color',clt(rel+1,:));
                        ylim([1, 4]) % rating range
                    end
                    xticks(abs_diff)
                end
            end

            %% model fitting

            % general setting for all models
            model.num_run        = num_run;
            model.num_sec         = 10; % number of samples in the parameter space, must be larger than num_run
            model.x               = (-512:1:512) * deg_per_px;
            model.sA              = sA;
            model.sV              = model.sA;
            model.num_rep         = num_rep;
            model.num_SD          = 5;
            model.numBins_A       = 15;
            model.numBins_V       = 15;
            model.modality        = {'A','V'};
            model.strategy_loc    = 'MA';

            OPTIONS.TolMesh = 1e-3;
            
            data.org_resp         = org_loc;
            data.org_conf         = org_conf;
            data.sigMotor         = fixP.sigMotor;
            saveData{d,i}         = data;

            currModel = str2func('nllBimodal');

            % switch confidence strategies
            model.strategy_conf         = ds_conf{d};

            % initiate
            model.mode                  = 'initiate';
            Val = currModel([], model, data);

            % optimize
            model.mode                  = 'optimize';
            NLL                         = NaN(1, model.num_run);
            estP                        = NaN(model.num_run, Val.num_para);

            % % test using ground truth parameters
            %             p = GT{d};
            %             test = currModel(p, model, data);

            parfor n              = 1:model.num_run

                tempModel             = model;
                tempVal               = Val;
                tempFunc              = currModel;

                [estP(n,:),NLL(n),~,~,traj(n)]    = bads(@(p) tempFunc(p, model, data),...
                    Val.init(n,:), Val.lb,...
                    Val.ub, Val.plb, Val.pub, [], OPTIONS);

                disp(estP(n,:))

            end

            % find the parameter with the least NLL
            [minNLL, best_idx]    = min(NLL);
            bestP                 = estP(best_idx, :);

            % save all fitting results
            saveModel{d,i}.estP = estP;
            saveModel{d,i}.NLL = NLL;
            saveModel{d,i}.bestP = bestP;
            saveModel{d,i}.minNLL = minNLL;

            % predict using the best-fitting parameter
            model.mode            = 'predict';
            tmpPred               = currModel(bestP, model, data);
            pred{d,i}             = tmpPred;

        end
    end

    save(flnm, 'saveModel','saveData', 'pred')

end