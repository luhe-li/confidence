
muP = 0;
sigP = 5;

% Define parameters for the first Gaussian distribution
mV = 0; % mean
sigV = 1; % standard deviation

% Define parameters for the second Gaussian distribution
mA = 5; % mean
sigA = 1.5; % standard deviation

% common cause
c = 0.65;

% Create the x range for the PDF
x = linspace(mV - 4*sigV, mA + 4*sigA, 1000);

% Calculate the PDFs
pdf_V = normpdf(x, mV, sigV);
pdf_V = pdf_V ./ sum(pdf_V);
pdf_A = normpdf(x, mA, sigA);
pdf_A = pdf_A ./ sum(pdf_A);

% Calculate reliability (inverse of variance)
reliabilityV = 1 / sigV^2;
reliabilityA = 1 / sigA^2;
reliabilityP = 1 / sigP^2;

% p(shat_a, c=2)
shat_A_C2_mean = (mA./sigA^2 + muP./sigP^2)./(1/sigA^2+1/sigP^2);
shat_A_C2_var = 1/(reliabilityP + reliabilityA + reliabilityV);
shat_A_C2_pdf = normpdf(x, shat_A_C2_mean, sqrt(shat_A_C2_var));
shat_A_C2_pdf = shat_A_C2_pdf ./ sum(shat_A_C2_pdf);

% p(shat_a, c=1)
shat_A_C1_mean = (reliabilityV * mV + reliabilityA * mA + reliabilityP * muP) / (reliabilityV + reliabilityA + reliabilityP);
shat_A_C1_variance = 1 / (reliabilityV + reliabilityA + reliabilityP);
shat_A_C1_pdf = normpdf(x, shat_A_C1_mean, sqrt(shat_A_C1_variance));
shat_A_C1_pdf = shat_A_C1_pdf ./ sum(shat_A_C1_pdf);

lw = 2;
opac = 0.5;

%%
% Plot the PDFs
figure;

subplot(3,1,1);
% Plot the grey shade under the first Gaussian curve
fill([x, fliplr(x)], [shat_A_C1_pdf, zeros(size(shat_A_C1_pdf))], ones(1,3).*opac, 'EdgeColor', 'none');
hold on;
% Plot the Gaussian curves
plot(x, pdf_V, '--','Color','#B4222D','LineWidth',lw);
[value,ind] = max(pdf_V);
plot([x(ind), x(ind)], [0, value], '--','Color','#B4222D','LineWidth',lw);

plot(x, pdf_A, '--','Color','#2270B0','LineWidth',lw);
[value,ind] = max(pdf_A);
plot([x(ind), x(ind)], [0, value], '--','Color','#2270B0','LineWidth',lw);

plot(x, shat_A_C1_pdf, '-','LineWidth',lw,'Color','#882574');
[value,ind] = max(shat_A_C1_pdf);
plot([x(ind), x(ind)], [0, value], '-', 'LineWidth', lw,'Color','#882574');
% Add labels and legend
% legend('','${m_V}$','${m_A}$','$\hat{s}$', 'Interpreter', 'latex');
hold off;
set(gca, 'XTick', [], 'XTickLabel', [],'YTick', [], 'YTickLabel', []);
xlim([min(x),max(x)]);

%%

subplot(3,1,2);
% Plot the grey shade under the first Gaussian curve
fill([x, fliplr(x)], [shat_A_C2_pdf, zeros(size(shat_A_C2_pdf))], ones(1,3).*opac, 'EdgeColor', 'none');
hold on;
lw = 2;

plot(x, pdf_A, '--','Color','#2270B0','LineWidth',lw);
[value,ind] = max(pdf_A);
plot([x(ind), x(ind)], [0, value], '-','Color','#2270B0','LineWidth',lw);

plot(x, shat_A_C2_pdf, 'Color','#2270B0','LineWidth',lw);
[value,ind] = max(shat_A_C2_pdf);
plot([x(ind), x(ind)], [0, value], 'Color','#2270B0', 'LineWidth', lw);

% Add labels and legend
% legend('','${m_V}$','${m_A}$','$\hat{s}$', 'Interpreter', 'latex');
hold off;
set(gca, 'XTick', [], 'XTickLabel', [],'YTick', [], 'YTickLabel', []);
xlim([min(x),max(x)]);

%%

subplot(3,1,3);
full_post = shat_A_C2_pdf .* (1-c) + shat_A_C1_pdf .* c;
full_post = full_post ./ sum(full_post);
% Plot the grey shade under the first Gaussian curve
fill([x, fliplr(x)], [full_post, zeros(size(full_post))], ones(1,3).*opac, 'EdgeColor', 'none');
hold on;
lw = 2;
% Plot the Gaussian curves
% 
% plot(x, shat_A_C2_pdf, '--','Color','#2270B0','LineWidth',lw);
[~,indA] = max(shat_A_C2_pdf);
% plot([x(indA), x(indA)], [0, value], '--','Color','#2270B0','LineWidth',lw);
% 
% plot(x, shat_A_C1_pdf, 'k--','LineWidth',lw);
[~,indC] = max(shat_A_C1_pdf);
% plot([x(indC), x(indC)], [0, value], 'k--', 'LineWidth', lw);

plot(x, full_post, '-','Color','k','LineWidth',lw);
[value,ind] = max(full_post);

MAestLoc = x(indA).* (1-c) + x(indC).* c;
plot([MAestLoc, MAestLoc], [0, value*1.1], '-', 'Color','k','LineWidth', lw);
% Add labels and legend
% legend('','${m_V}$','${m_A}$','$\hat{s}$', 'Interpreter', 'latex');
hold off;
set(gca, 'XTick', [], 'XTickLabel', [],'YTick', [], 'YTickLabel', []);
xlim([min(x),max(x)]);

out_dir               = fullfile(pwd, mfilename);
if ~exist(out_dir,'dir'); mkdir(out_dir); end
saveas(gca,fullfile(out_dir,'cau_inf_model'),'pdf')






