%% localization responses as a function of stimulus location 
idx_loc_lb      = arrayfun(@(idx) find(~isnan(pC1_resp_lenD(1,1,:,idx)),1),1:lenD);
                                                      %[4,3,2,1,1,1,1];
idx_loc_ub      = arrayfun(@(idx) ceil(lenD/2)-find(~isnan(...
                      flipud(squeeze(pC1_resp_lenD(1,1,:,idx)))),1), 1:lenD)+1;
                                                      %[4,4,4,4,3,2,1];             
idx_sloc_lb     = fliplr(idx_loc_lb(1,:));            %[1,1,1,1,2,3,4];
idx_sloc_ub     = fliplr(idx_loc_ub(1,:));            %[1,2,3,4,4,4,4];
ybd_locR        = [min(mean_LocResp(:))-max(SD_LocResp(:)),...
                    max(mean_LocResp(:))+max(SD_LocResp(:))]; 
y_ticks_locR    = V_loc;
lgd_label       = cell(1,lenD);for l=1:lenD;lgd_label{l}=num2str(AV_discrepancy(l));end
x_lbl           = {'Auditory','Visual'};
cMap_locR(1,:,:)= [0.2706, 0.3059, 0.5529;0.3216, 0.4824, 0.6471;0.3725, 0.6549, 0.7412;...
                   0.0588, 0.7373, 0.7451;0.3529, 0.7373, 0.4588;0.5686, 0.8039, 0;...
                   0.7373, 0.8471, 0];
cMap_locR(2,:,:)= [0.4275, 0.1059, 0.1647; 0.6157, 0.1490, 0.3843; 0.7216, 0.4078, 0.4824;...
                   0.8667, 0.5686, 0.5725; 0.9451, 0.6353, 0.5333; 0.9490, 0.7098, 0.5020;...
                   0.9608, 0.7922, 0.4431];
subplot_idx     = [1,3;2,4];
                
if bool_plt(2) == 1
    for m = 1:lenM
        figure(5+m)
        %x boundaries change based on modality
        xbd     = [eval([modality{m}, '_loc(1)'])-3, eval([modality{m},'_loc(end)'])+3]; 
        x_ticks = round(eval([modality{m}, '_loc']),1);
        for i = 1:lenC
            for j = 1:lenP
                subplot(lenC, lenP, subplot_idx(i,j))
                addBackground(xbd, ybd_locR, [xbd(1),x_ticks, xbd(end)], ...
                    [ybd_locR(1),y_ticks_locR,ybd_locR(end)]);
                plot(eval([modality{m},'_loc']), V_loc,'Color',ones(1,3).*0.5,...
                    'lineWidth',lw,'lineStyle','--'); hold on
                for k = 1:lenD
                    Eb(k) = errorbar(eval([modality{m},'_loc(idx_sloc_lb(k):idx_sloc_ub(k))']), ...
                        squeeze(mean_LocResp(i,j,m,idx_loc_lb(k):idx_loc_ub(k),k)),...
                        squeeze(SD_LocResp(i,j,m,idx_loc_lb(k):idx_loc_ub(k),k)),'-s',...
                        'Color',cMap_locR(m,k,:),'lineWidth', lw,'MarkerSize',15,...
                        'MarkerEdgeColor',cMap_locR(m,k,:),'MarkerFaceColor',...
                        cMap_locR(m,k,:)); hold on
                end
                hold off; box off; axis square;
                xticks(x_ticks); xlim(xbd); xlabel([x_lbl{m},' stimulus location (deg)']);
                yticks(y_ticks_locR); ylim(ybd_locR); 
                ylabel(sprintf([x_lbl{m},' mean localization\nresponses (deg)']));
                if i==1 && j == 2
                    lgd = legend(Eb,lgd_label,'Location','northeast','FontSize',fs_lgds,...
                        'Orientation','horizontal');legend boxoff; htitle = get(lgd,'Title'); 
                    set(htitle,'String',['Spatial discrepancy (V-A, deg)']);
                end
                if i==1 && j==1; text(xbd(1)+0.2,ybd_locR(end)-1,subjI,'FontSize',fs_lgds); end
                %title(['Condition: ', cond{i}, ', phase: ', phase{j}],'FontSize',fs_lgds);
                set(gca,'FontSize',fs_lbls);
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.85]);
            end
        end
        set(gcf,'PaperUnits','centimeters','PaperSize',[40 40]);
        %saveas(gcf, ['mean_LocResp',modality{m},'_',cond{i}, '_',subjI], 'pdf'); 
    end
end
