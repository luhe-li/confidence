%bimodal A
figure; imagesc(r_bimodal_param_btst); colormap(winter);colorbar;
xticks(2:2:10); yticks(200:200:1000); 
caxis([min(r_bimodal_param_btst(:)), max(r_bimodal_param_btst(:))]);
xlabel('Simulated response trial'); ylabel('Parametric bootstrapping run');
title('$r_{AV_l,A}$','interpreter','latex');set(gca,'FontSize',20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.25, 0.7]); 
set(gcf,'PaperUnits','centimeters','PaperSize',[20 35]); 
%saveas(gcf, 'SimulatedResp_A_paramBtst','pdf');

%unimodal A r_uni_param_btst
figure; imagesc(r_uni_param_btst_kj); colormap(winter);colorbar;
xticks(5:5:30); yticks(200:200:1000); 
caxis([min(r_bimodal_param_btst(:)), max(r_bimodal_param_btst(:))]);
xlabel('Simulated response trial'); ylabel('Parametric bootstrapping run');
title('$r_{A}$','interpreter','latex');set(gca,'FontSize',20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.7]); 
set(gcf,'PaperUnits','centimeters','PaperSize',[40 35]); 
%saveas(gcf, 'SimulatedResp_uniA_paramBtst','pdf');

figure; imagesc(squeeze(r_uni_param_btst(2,1,:))); colormap(winter);colorbar;
yticks(200:200:1000); xticks([]);
caxis([min(r_bimodal_param_btst(:)), max(r_bimodal_param_btst(:))]);
ylabel('Parametric bootstrapping run');
title('$\overline{r}_{A}$','interpreter','latex');set(gca,'FontSize',20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.12, 0.7]); 
set(gcf,'PaperUnits','centimeters','PaperSize',[10 35]); 
%saveas(gcf, 'SimulatedResp_uniA_mean_paramBtst','pdf');

%diff
figure;imagesc(VE_param_btst_nanmean');colorbar
xticks(1:1:7); xticklabels(-24:8:24); yticks(200:200:1000); 
caxis([min(VE_param_btst_nanmean(:)), max(VE_param_btst_nanmean(:))]);
xlabel('Spatial discrepancy (deg)'); ylabel('Parametric bootstrapping run');
set(gca,'FontSize',20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.7]); 
set(gcf,'PaperUnits','centimeters','PaperSize',[25 35]); 
saveas(gcf, 'A_VE_paramBtst','pdf');

figure;imagesc(repmat([-1,-1,-1,1,1,1,1],[1000,1]).*VE_param_btst_nanmean');colorbar
xticks(1:1:7); xticklabels(-24:8:24); yticks(200:200:1000); 
caxis([min(VE_param_btst_nanmean(:)), max(VE_param_btst_nanmean(:))]);
xlabel('Spatial discrepancy (deg)'); ylabel('Parametric bootstrapping run');
set(gca,'FontSize',20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.7]); 
set(gcf,'PaperUnits','centimeters','PaperSize',[25 35]); 
saveas(gcf, 'A_VE_corr_paramBtst','pdf');

figure;imagesc(VE_abs_param_btst_sort');colorbar
xticks(1:1:4); xticklabels(0:8:24); yticks(200:200:1000); 
caxis([min(VE_param_btst_nanmean(:)), max(VE_param_btst_nanmean(:))]);
xlabel('Absolute spatial discrepancy (deg)'); ylabel('Parametric bootstrapping run');
set(gca,'FontSize',20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.25, 0.7]); 
set(gcf,'PaperUnits','centimeters','PaperSize',[20 35]); 
saveas(gcf, 'A_VE_abs_paramBtst','pdf');

