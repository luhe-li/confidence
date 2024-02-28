clear; clc;
sub_slc = 1;
ses_slc = 'V';
total_num_rep = 20;
exp = 'uniLoc';

% manage path
cur_dir      = pwd;
[project_dir, ~]= fileparts(fileparts(cur_dir));
% out_dir      = fullfile(cur_dir, 's1Fig');
data_dir     = fullfile(project_dir,'data',exp);
% if ~exist(out_dir,'dir') mkdir(out_dir); end

% organize data
flnm        = sprintf('uniLoc_sub%i_ses-V', sub_slc);

load(fullfile(data_dir, flnm))
%%
% sigVperTrial = NaN(1,length(seq));
confTrial = NaN(1,length(seq));
errTrial = NaN(1,length(seq));

for i = 1:length(seq)
% sigVperTrial(i) = std(Resp(i).vStimDotsCoor(1,:),[],2);
confTrial = org_conf(:);
errTrial(i) = Resp(i).response_pixel - 512 - Resp(i).target_pixel;
end

%%
figure
h1 = histogram(abs(errTrial(~~confTrial)),50);
hold on
h2 = histogram(abs(errTrial(~confTrial)),50);
hold off
legend('Confident','Not Confident')
ylabel('Count')
xlabel('Absolute error (pixel)')
%%
figure
h1 = histogram(sigVperTrial(~~confTrial),50);
hold on
h2 = histogram(sigVperTrial(~confTrial),50);
hold off
legend('Confident','Not Confident')
ylabel('Count')
xlabel('\sigma_{v} (pixel)')