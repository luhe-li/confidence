function plot_spread_VE(bi_resp,aud_locs,raw_diff,remapped_vis_locs,cue,reliability)
lw = 3;
cueLabel = [{'Audio'},{'Visual'}];
resp_rel = squeeze(mean(bi_resp(:,:,cue,reliability,:),5));
resp_rel_std = std(squeeze(bi_resp(:,:,cue,reliability,:)),[],3);

% uni resp unneeded. if use, add an argument in function input
% if size(uni_resp,2) == 6
%     a_uni = squeeze(mean(uni_resp(1,2:5,:),3));
% elseif size(uni_resp,2) == 4
%     a_uni = squeeze(mean(uni_resp(1,:,:),3));
% end
% a_uni_mat = repmat(a_uni,[4,1])';
switch cue
    case 1
        stim_locs_mat = repmat(aud_locs,[4,1])';
        index = aud_locs - aud_locs';
        diff_vec = unique(round(index));
        plot(aud_locs, remapped_vis_locs,'Color',ones(1,3).*0.5,...
                    'lineWidth',lw,'lineStyle','--'); hold on
    case 2
        stim_locs_mat = repmat(remapped_vis_locs,4,1);
        index = remapped_vis_locs - remapped_vis_locs';
        diff_vec = unique(round(index));
        plot(remapped_vis_locs, remapped_vis_locs,'Color',ones(1,3).*0.5,...
                    'lineWidth',lw,'lineStyle','--'); hold on
end


cMap_locR(1,:,:)= [0.2706, 0.3059, 0.5529; 0.30885, 0.438275, 0.62355; 0.34705, 0.56865, 0.69415; 0.294075, 0.6755, 0.742175; 0.0588, 0.7373, 0.7451; 0.279375, 0.7373, 0.530375; 0.46075, 0.7706, 0.2294; 0.610775, 0.8147, 0; 0.7373, 0.8471, 0];
cMap_locR(2,:,:)= [0.4275, 0.1059, 0.1647; 0.56865, 0.138225, 0.3294; 0.66865, 0.2784, 0.43335; 0.757875, 0.448, 0.504925; 0.8667, 0.5686, 0.5725; 0.9255, 0.618625, 0.5431; 0.94705, 0.67255, 0.51765; 0.95195, 0.7304, 0.487275; 0.9608, 0.7922, 0.4431];

xline(0,'--');
yline(0,'--');
plot(-15:15, -15:15,'k--')
for i = 1:length(diff_vec)
    diff_i = diff_vec(i);
    Eb(i) = errorbar(stim_locs_mat(round(index) == diff_i),resp_rel(round(index) == diff_i),resp_rel_std(round(index) == diff_i),'-o' ...
        ,'Color',cMap_locR(cue,i,:),'lineWidth', lw,'MarkerSize',15,...
         'MarkerEdgeColor',cMap_locR(cue,i,:),'MarkerFaceColor',...
         cMap_locR(cue,i,:));
    hold on
end
hold off
xlabel(sprintf([cueLabel{cue} ' stimulus location (deg)']),'FontSize',20);
ylabel(sprintf([cueLabel{cue} '  mean localization\nresponses (deg)']),'FontSize',20);
lgd_label       = cell(1,length(raw_diff));for l=1:length(raw_diff);lgd_label{l}=num2str(round(raw_diff(l)));end
lgd = legend(Eb,lgd_label,'Location','northoutside','FontSize',15,...
    'Orientation','horizontal');legend boxoff; htitle = get(lgd,'Title');
set(htitle,'String','Spatial discrepancy (V-A, deg)');
set(gca,'FontSize',20);
end