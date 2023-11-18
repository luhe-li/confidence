clear all; close all; clc

%% different combinations of variability of sensory measurements
sA      = -2;
sV      = 2;
nSteps  = 100;
sigma_A = linspace(0.1, 5,nSteps);
sigma_V = linspace(0.01,2, nSteps);
mu_P    = 0;
sigma_P = 100;
pCommon = 0.37;
nSims   = 1e5;
post_C1 = NaN(nSteps, nSims);

for i = 1:nSteps
    %constants for fitting the unity judgment
    CI.J_A            = sigma_A(i)^2;
    CI.J_V            = sigma_V(i)^2;
    CI.J_P            = sigma_P^2;
    CI.constC1        = CI.J_A*CI.J_V + CI.J_A*CI.J_P + CI.J_V*CI.J_P;
    CI.constC1_shat   = 1/CI.J_A  + 1/CI.J_V + 1/CI.J_P;
    CI.constC2_1      = CI.J_A + CI.J_P;
    CI.constC2_1_shat = 1/CI.J_A + 1/CI.J_P;
    CI.constC2_2      = CI.J_V + CI.J_P;
    CI.constC2_2_shat = 1/CI.J_V + 1/CI.J_P; 
    
    mA = randn(1, nSims).*sigma_A(i) + sA;
    mV = randn(1, nSims).*sigma_V(i) + sV;
    [post_C1(i,:), ~, ~, ~] = calculatePostC1C2(mA, mV, mu_P, CI, pCommon);
end

mean_post_C1 = mean(post_C1,2);
post_C1_sort = sort(post_C1,2);
post_C1_lb   = post_C1_sort(:, round(nSims*0.16));
post_C1_ub   = post_C1_sort(:, round(nSims*0.84));
plotPostC1('sigma', sigma_A, sigma_V, sA, sV, post_C1_lb, ...
    post_C1_ub, mean_post_C1, nSteps, 1);

%% the variabilities remain unchanged, but the spatial offset between the 
%auditory and visual measurements changes
sA      = -5;
sV      = linspace(-5, 8, nSteps);
sigma_A = 2.5;
sigma_V = 0.5;
post_C1 = NaN(nSteps, nSims);

%constants for fitting the unity judgment
CI.J_A            = sigma_A^2;
CI.J_V            = sigma_V^2;
CI.J_P            = sigma_P^2;
CI.constC1        = CI.J_A*CI.J_V + CI.J_A*CI.J_P + CI.J_V*CI.J_P;
CI.constC1_shat   = 1/CI.J_A  + 1/CI.J_V + 1/CI.J_P;
CI.constC2_1      = CI.J_A + CI.J_P;
CI.constC2_1_shat = 1/CI.J_A + 1/CI.J_P;
CI.constC2_2      = CI.J_V + CI.J_P;
CI.constC2_2_shat = 1/CI.J_V + 1/CI.J_P; 
for i = 1:nSteps    
    mA = randn(1, nSims).*sigma_A + sA;
    mV = randn(1, nSims).*sigma_V + sV(i);
    [post_C1(i,:), ~, ~, ~] = calculatePostC1C2(mA, mV, mu_P, CI, pCommon);
end

mean_post_C1 = mean(post_C1,2);
post_C1_sort = sort(post_C1,2);
post_C1_lb   = post_C1_sort(:, round(nSims*0.16));
post_C1_ub   = post_C1_sort(:, round(nSims*0.84));
plotPostC1('distance', sigma_A, sigma_V, sA, sV, post_C1_lb, ...
    post_C1_ub, mean_post_C1, nSteps, 1);

%%
function [Post_C1, Post_C2, L_C1, L_C2] = calculatePostC1C2(X1, X2, mu_P, CI, pC1)
    %likelihood of a common cause and seperate causes
    L_C1     = 1/(2*pi*sqrt(CI.constC1))*exp(-0.5*((X1 - X2).^2.*CI.J_P +...
                (X1 - mu_P).^2*CI.J_V + (X2 - mu_P).^2.*CI.J_A)./CI.constC1);
    L_C2     = 1/(2*pi*sqrt(CI.constC2_1*CI.constC2_2))*exp(-0.5*...
                ((X1 - mu_P).^2./CI.constC2_1+(X2 - mu_P).^2./CI.constC2_2)); 
    normTerm = L_C1.*pC1 + L_C2.*(1-pC1);
    %posterior of a common cause
    Post_C1  = L_C1.*pC1./normTerm;
    Post_C2  = 1 - Post_C1;
end

function plotPostC1(var_slc, sigma_A, sigma_V, sA, sV, post_C1_lb, ...
    post_C1_ub, mean_post_C1, nSteps, bool_save)
    slc_idx = [20,40,80];
    xx      = linspace(-10,10,1e3);
    cb      = [65,105,225]./255; 
    cr      = [0.85, 0.23, 0.28];
    cmap    = [255, 188,107; 113,159,94; 100, 85,139]./255;


    figure
    for i = 1:length(slc_idx)
        subplot(2,3,i)
        patch([xx(1), xx(end), xx(end), xx(1)], [0, 0, 1, 1], cmap(i,:),...
            'FaceAlpha', 0.1, 'EdgeColor', cmap(i,:),'lineWidth',2); hold on
        if strcmp(var_slc, 'sigma')
            plot(xx, normpdf(xx, sA, sigma_A(slc_idx(i))), 'Color',cb,'lineWidth',2); 
            plot(xx, normpdf(xx, sV, sigma_V(slc_idx(i))), 'Color',cr,'lineWidth',2); 
        elseif strcmp(var_slc, 'distance')
            plot(xx, normpdf(xx, sA, sigma_A), 'Color',cb,'lineWidth',2); 
            plot(xx, normpdf(xx, sV(slc_idx(i)), sigma_V), 'Color',cr,'lineWidth',2);
        end
        ylim([0, 1]); xlim([xx(1), xx(end)]); box off; xticks([]); yticks([]);
    end

    subplot(2,3,[4,5,6])
    patch([1:nSteps, nSteps:-1:1], [post_C1_lb', flipud(post_C1_ub)'],...
        'k','FaceAlpha',0.1,'EdgeColor','none'); hold on
    plot(1:nSteps, mean_post_C1,'k','lineWidth',3); 
    for i = 1:length(slc_idx)
        scatter(slc_idx(i), mean_post_C1(slc_idx(i)), 100, ...
            'MarkerFaceColor',cmap(i,:),'Marker','d','markerEdgeColor','k'); 
    end
    hold off; box off; xticks([1, 50, nSteps]); yticks(0:0.5:1); xlim([1, nSteps]);
    if strcmp(var_slc, 'sigma')
        xticklabels({'Low','Medium','High'});
        xlabel('Variability of sensory measurements'); 
    elseif strcmp(var_slc, 'distance')
        xticklabels({'None','Close','Far'});
        xlabel('Distance between measurements');  
    end
    ylabel(sprintf('Posterior probability \nof common cause'));
    set(gca,'FontSize',15);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.35, 0.45]);
    set(gcf,'PaperUnits','centimeters','PaperSize',[25 25]);
    if bool_save
        if strcmp(var_slc, 'sigma'); saveas(gcf, 'Fig6A', 'pdf'); 
        elseif strcmp(var_slc, 'distance'); saveas(gcf, 'Fig6B', 'pdf'); 
        end
    end
end
