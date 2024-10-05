
clear; clc;
model_slc = 5;

%% model info

rng('Shuffle');
specifications = {'MA, full posterior, global','MA, full posterior, local','MA, gaussian posterior','MS','PM'};
folders = {'MA_optimal', 'MA_local', 'MA_gauss', 'MS','PM'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% manage path

restoredefaultpath;
data_dir = (fullfile(pwd, 'param_recovery'));
out_dir = fullfile(pwd, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load and clean data

for mm = 1:numel(model_slc)

    i_model = model_slc(mm);
    file_pattern = sprintf('%s-*.mat', folders{i_model});
    files = dir(fullfile(data_dir, file_pattern));
    if ~isempty(files)
        flnm = files(1).name;
        load(fullfile(data_dir, flnm));
        disp(['Loaded file: ' flnm]);
    else
        disp('No files matching the pattern found.');
    end

    % use GT from 1 sample because it's the same
    GT(mm, :) = sim_data(model_slc,1,1).gt;

    for ss = 1:100
        
        best_p(mm, ss, :) = fits(model_slc,1,ss).best_p;

    end

end

%% plot
for mm = 1:numel(model_slc)
    figure;
    hold on;
    n_para = size(GT, 2);
    
    for pp = 1:n_para
        subplot(1, n_para, pp);
        histogram(squeeze(best_p(mm, :, pp)), 'Normalization', 'pdf');
        hold on;
        xline(GT(mm, pp), 'r', 'LineWidth', 2);
        title(['Parameter ' num2str(pp)]);
        xlabel('Ground-truth');
        ylabel('Count');
    end
    
    sgtitle(['Model ' num2str(mm)]);
end

