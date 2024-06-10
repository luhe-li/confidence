clear; %close all; 
rng('shuffle');

sub_slc = 13;
ses_slc = 1:2; % bimodal sessions

models = {'Heuristic','Suboptimal','Optimal'};

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(fileparts(cur_dir));
out_dir               = fullfile(cur_dir, mfilename);
result_dir            = fullfile(cur_dir, 'modelFit');
if ~exist(out_dir,'dir') mkdir(out_dir); end
addpath(genpath(fullfile(project_dir, 'func')))
addpath(genpath(result_dir))

%% plot set up

lw = 1;
pred_lw = 2;
fontSZ = 15;
titleSZ = 20;
dotSZ = 80;
clt = [30, 120, 180; % blue
    227, 27, 27;  % dark red
    repmat(125, 1, 3)]./255;

%% predict data

% load best fitting parameters
flnm = sprintf('fitResults_sub%i_ses%i-%i.mat', sub_slc, min(ses_slc), max(ses_slc));
files = dir(fullfile(result_dir, flnm));
load(files(end).name);

% winning model?
[~, d] = min([saveConfModel{1}.minNLL, saveConfModel{2}.minNLL, saveConfModel{3}.minNLL]);
d = 3;

% use the best-fitting parameter and winning model
p = saveConfModel{d}.bestP;
model.mode = 'predict';
model.model_slc             = d;
model.strategy_conf         = models{d};

% increase trial number for prediction
model.uni_nrep = 1e3;
model.bi_nrep = 1e3;
pred = nllUniBiConf(p, model, data);
disp(p)

%% plot unimodal estvar

uni_est = reshape(pred.uni_est_var, 3,[]);
figure; hold on
histogram(uni_est(1,:))
histogram(uni_est(2,:))
histogram(uni_est(3,:))
legend({'A','V1','V2'})
cA = [p(2),p(2)+p(3),p(2)+p(3)+p(4)]; cV = [p(5),p(5)+p(6),p(5)+p(6)+p(7)];

for c = 1:3

    xline(cA(c));
    xline(cV(c),'--');
end


%% plot unimodal estvar

uni_est = reshape(pred.uni_est_var, 3,[]);
figure; hold on
histogram(uni_est(1,:))
histogram(uni_est(2,:))
histogram(uni_est(3,:))
legend({'A','V1','V2'})
cA = [p(2),p(2)+p(3),p(2)+p(3)+p(4)]; cV = [p(5),p(5)+p(6),p(5)+p(6)+p(7)];

for c = 1:3

    xline(cA(c));
    xline(cV(c),'--');

end

%% plot bimodal estvar

sA = [-10, -2, 2, 10];
diff =  sA - sA';
udiff = unique(diff);
udiff = udiff(5:end);
pair_slc = [2,2;3,2;2,1;3,1;4,1];
npair = size(pair_slc,1);

figure
sgtitle('A')
set(gcf,'Position',[0,0,2400,300])

for pp = 1:npair

    subplot(1,npair, pp);
    hold on
    % xlim([0, 10])
    title(sprintf('discrepancy: %i', udiff(pp)))

    aa = pair_slc(pp,1);
    vv = pair_slc(pp,2);

    for  rel = 1:2

       bi_est = squeeze(pred.bi_est_var(aa,vv,1,rel,:));
       histogram(bi_est,'NumBins',25)

    end

    % for c = 1:3

        % xline(cA(c));
        % xline(cV(c),'--');
    % end

end


%% plot bimodal conf resp

sA = [-10, -2, 2, 10];
diff =  sA - sA';
udiff = unique(diff);
udiff = udiff(5:end);
pair_slc = [2,2;3,2;2,1;3,1;4,1];
npair = size(pair_slc,1);

figure
sgtitle('A')
set(gcf,'Position',[0,0,2400,300])

for pp = 1:npair

    subplot(1,npair, pp);
    hold on
    % xlim([0, 10])
    title(sprintf('discrepancy: %i', udiff(pp)))

    aa = pair_slc(pp,1);
    vv = pair_slc(pp,2);

    for  rel = 1:2

       bi_conf = squeeze(pred.bi_conf(aa,vv,1,rel,:));
       histogram(bi_conf(:))

    end
    xticks(0:5)
    xlim([0,5])

end

%%
figure
sgtitle('V')
set(gcf,'Position',[0,0,2400,300])

for pp = 1:npair

    subplot(1,npair, pp);
    hold on
    % xlim([0, 10])
    title(sprintf('discrepancy: %i', udiff(pp)))

    aa = pair_slc(pp,1);
    vv = pair_slc(pp,2);

    for  rel = 1:2

       bi_est = squeeze(pred.bi_est_var(aa,vv,2,rel,:));
       histogram(bi_est,'NumBins',25)

    end

    for c = 1:3

        % xline(cA(c));
        % xline(cV(c),'--');
    end

end
