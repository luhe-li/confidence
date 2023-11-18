%set breakpoint at plotModelFits.m/line 204
figure; imagesc(rBinary_param_btst); colormap(hot); colorbar
xticks(10:10:40); yticks(200:200:1000); caxis([0,1]);
xlabel('Simulated response trial'); ylabel('Parametric bootstrapping run');
title('$I_{C=1}$','interpreter','latex');set(gca,'FontSize',20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.7]); 
set(gcf,'PaperUnits','centimeters','PaperSize',[30 35]); 
saveas(gcf, 'SimulatedResp_unity_paramBtst','pdf');

figure; imagesc(propC1_param_btst_ij'); colormap(hot); colorbar
xticks(1:1:7); xticklabels(-24:8:24); yticks(200:200:1000); caxis([0,1]);
xlabel('Spatial discrepancy (deg)'); ylabel('Parametric bootstrapping run');
set(gca,'FontSize',20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.7]); 
set(gcf,'PaperUnits','centimeters','PaperSize',[25 35]); 
saveas(gcf, 'SimulatedProportionC1_unity_paramBtst','pdf');

figure; imagesc(propC1_abs_param_btst'); colormap(hot); colorbar
xticks(1:1:4); xticklabels(0:8:24); yticks(200:200:1000); caxis([0,1]);
xlabel('Absolute spatial discrepancy (deg)'); ylabel('Parametric bootstrapping run');
set(gca,'FontSize',20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.25, 0.7]); 
set(gcf,'PaperUnits','centimeters','PaperSize',[20 35]); 
saveas(gcf, 'SimulatedProportionC1_unity_abs_paramBtst','pdf');

figure; imagesc(propC1_abs_param_btst_sort'); colormap(hot); colorbar
xticks(1:1:4); xticklabels(0:8:24); yticks(200:200:1000); caxis([0,1]);
xlabel('Absolute spatial discrepancy (deg)'); ylabel('Parametric bootstrapping run');
set(gca,'FontSize',20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.25, 0.7]); 
set(gcf,'PaperUnits','centimeters','PaperSize',[20 35]); 
saveas(gcf, 'SimulatedProportionC1_unity_abs_sorted_paramBtst','pdf');

