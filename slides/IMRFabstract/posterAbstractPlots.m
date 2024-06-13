% Define parameters for the first Gaussian distribution
mu1 = 0; % mean
sigma1 = 1; % standard deviation

% Define parameters for the second Gaussian distribution
mu2 = 5; % mean
sigma2 = 1.5; % standard deviation

% common cause
c = 0.65;
% Create the x range for the PDF
x = linspace(mu1 - 4*sigma1, mu2 + 4*sigma2, 1000);

% Calculate the PDFs
pdf_V = normpdf(x, mu1, sigma1);
pdf_V = pdf_V ./ sum(pdf_V);
pdf_A = normpdf(x, mu2, sigma2);
pdf_A = pdf_A ./ sum(pdf_A);
% Calculate reliability (inverse of variance)
reliability1 = 1 / sigma1^2;
reliability2 = 1 / sigma2^2;

% Weighted combination of the two PDFs
combined_mean = (reliability1 * mu1 + reliability2 * mu2) / (reliability1 + reliability2);
combined_variance = 1 / (reliability1 + reliability2);
combined_pdf = normpdf(x, combined_mean, sqrt(combined_variance));
combined_pdf = combined_pdf ./ sum(combined_pdf);
lw = 2;
opac = 0.5;
%%
% Plot the PDFs
figure;

subplot(3,1,1);
% Plot the grey shade under the first Gaussian curve
fill([x, fliplr(x)], [combined_pdf, zeros(size(combined_pdf))], ones(1,3).*opac, 'EdgeColor', 'none');
hold on;
% Plot the Gaussian curves
plot(x, pdf_V, '--','Color','#B4222D','LineWidth',lw);
[value,ind] = max(pdf_V);
plot([x(ind), x(ind)], [0, value], '--','Color','#B4222D','LineWidth',lw);

plot(x, pdf_A, '--','Color','#2270B0','LineWidth',lw);
[value,ind] = max(pdf_A);
plot([x(ind), x(ind)], [0, value], '--','Color','#2270B0','LineWidth',lw);

plot(x, combined_pdf, 'k-','LineWidth',lw);
[value,ind] = max(combined_pdf);
plot([x(ind), x(ind)], [0, value], 'k-', 'LineWidth', lw);
% Add labels and legend
% legend('','${m_V}$','${m_A}$','$\hat{s}$', 'Interpreter', 'latex');
hold off;
set(gca, 'XTick', [], 'XTickLabel', [],'YTick', [], 'YTickLabel', []);
xlim([min(x),max(x)]);




subplot(3,1,2);
% Plot the grey shade under the first Gaussian curve
fill([x, fliplr(x)], [pdf_A, zeros(size(pdf_A))], ones(1,3).*opac, 'EdgeColor', 'none');
hold on;
lw = 2;
% Plot the Gaussian curves
plot(x, pdf_V, '--','Color','#B4222D','LineWidth',lw);
[value,ind] = max(pdf_V);
plot([x(ind), x(ind)], [0, value], '--','Color','#B4222D','LineWidth',lw);

plot(x, pdf_A, '-','Color','#2270B0','LineWidth',lw);
[value,ind] = max(pdf_A);
plot([x(ind), x(ind)], [0, value], '-','Color','#2270B0','LineWidth',lw);

plot(x, combined_pdf, 'k--','LineWidth',lw);
[value,ind] = max(combined_pdf);
plot([x(ind), x(ind)], [0, value], 'k--', 'LineWidth', lw);

% Add labels and legend
% legend('','${m_V}$','${m_A}$','$\hat{s}$', 'Interpreter', 'latex');
hold off;
set(gca, 'XTick', [], 'XTickLabel', [],'YTick', [], 'YTickLabel', []);
xlim([min(x),max(x)]);





subplot(3,1,3);
full_post = pdf_A .* (1-c) + combined_pdf .* c;
full_post = full_post ./ sum(full_post);
% Plot the grey shade under the first Gaussian curve
fill([x, fliplr(x)], [full_post, zeros(size(full_post))], ones(1,3).*opac, 'EdgeColor', 'none');
hold on;
lw = 2;
% Plot the Gaussian curves

plot(x, pdf_A, '--','Color','#2270B0','LineWidth',lw);
[value,indA] = max(pdf_A);
plot([x(indA), x(indA)], [0, value], '--','Color','#2270B0','LineWidth',lw);

plot(x, combined_pdf, 'k--','LineWidth',lw);
[value,indC] = max(combined_pdf);
plot([x(indC), x(indC)], [0, value], 'k--', 'LineWidth', lw);

plot(x, full_post, '-','Color','#882574','LineWidth',lw);
[value,ind] = max(full_post);

MAestLoc = x(indA).* (1-c) + x(indC).* c;
plot([MAestLoc, MAestLoc], [0, value*1.1], '-', 'Color','#882574','LineWidth', lw);
% Add labels and legend
% legend('','${m_V}$','${m_A}$','$\hat{s}$', 'Interpreter', 'latex');
hold off;
set(gca, 'XTick', [], 'XTickLabel', [],'YTick', [], 'YTickLabel', []);
xlim([min(x),max(x)]);







