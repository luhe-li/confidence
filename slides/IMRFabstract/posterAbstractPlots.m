% Define parameters for the first Gaussian distribution
mu1 = 0; % mean
sigma1 = 1; % standard deviation

% Define parameters for the second Gaussian distribution
mu2 = 5; % mean
sigma2 = 1.5; % standard deviation

% Create the x range for the PDF
x = linspace(mu1 - 4*sigma1, mu2 + 4*sigma2, 1000);

% Calculate the PDFs and normalize
pdfV = normpdf(x, mu1, sigma1);
pdfV = pdfV ./ sum(pdfV);

pdfA = normpdf(x, mu2, sigma2);
pdfA = pdfA ./ sum(pdfA);

% Calculate reliability (inverse of variance)
reliability1 = 1 / sigma1^2;
reliability2 = 1 / sigma2^2;

% Weighted combination of the two PDFs
combined_mean = (reliability1 * mu1 + reliability2 * mu2) / (reliability1 + reliability2);
combined_variance = 1 / (reliability1 + reliability2);
combined_pdf = normpdf(x, combined_mean, sqrt(combined_variance));
combined_pdf = combined_pdf ./ sum(combined_pdf);

% Combination coefficient
cc = 0.5;
lw = 2;
opac = 0.5;
% Plot the PDFs

figure;
% Subplot 1
subplot(3,1,1);
fill([x, fliplr(x)], [combined_pdf, zeros(size(combined_pdf))],  ones(1,3) .* opac, 'EdgeColor', 'none');
hold on;

plot(x, pdfV, '--','Color','#B4222D','LineWidth',lw);
plot(x, pdfA, '--','Color','#2270B0','LineWidth',lw);
plot(x, combined_pdf, 'k-','LineWidth',lw);

% Add vertical lines
plot([mu1, mu1], [0, max(pdfV)], '--','Color','#B4222D','LineWidth',lw);
plot([mu2, mu2], [0, max(pdfA)], '--','Color','#2270B0','LineWidth',lw);
plot([combined_mean, combined_mean], [0, max(combined_pdf)], 'k-','LineWidth',lw);
xlim([-2.7,9]);
set(gca, 'XTick', [], 'XTickLabel', [],'YTick', [], 'YTickLabel', []);
hold off;

% Subplot 2
subplot(3,1,2);
fill([x, fliplr(x)], [pdfA, zeros(size(pdfA))],  ones(1,3) .* opac, 'EdgeColor', 'none');
hold on;
plot(x, pdfV, '--','Color','#B4222D','LineWidth',lw);
plot(x, pdfA, '-','Color','#2270B0','LineWidth',lw);
plot(x, combined_pdf, 'k--','LineWidth',lw);

% Add vertical lines
plot([mu1, mu1], [0, max(pdfV)], '--','Color','#B4222D','LineWidth',lw);
plot([mu2, mu2], [0, max(pdfA)], '-','Color','#2270B0','LineWidth',lw);
plot([combined_mean, combined_mean], [0, max(combined_pdf)], 'k--','LineWidth',lw);
xlim([-2.7,9]);

set(gca, 'XTick', [], 'XTickLabel', [],'YTick', [], 'YTickLabel', []);
hold off;

% Subplot 3
subplot(3,1,3);
post = pdfA .* (1-cc) + combined_pdf .* cc;
post = post ./ sum(post);
MAloc = mu2 .* (1-cc) + combined_mean .* cc;

fill([x, fliplr(x)], [post, zeros(size(post))], ones(1,3) .* opac, 'EdgeColor', 'none');
hold on;
plot(x, post, '-','Color','#702963','LineWidth',lw);
plot(x, pdfA, '--','Color','#2270B0','LineWidth',lw);
plot(x, combined_pdf, 'k--','LineWidth',lw);

plot([mu2, mu2], [0, max(pdfA)], '--','Color','#2270B0','LineWidth',lw);
plot([combined_mean, combined_mean], [0, max(combined_pdf)], 'k--','LineWidth',lw);
plot([MAloc, MAloc], [0, max(post).* 1.5], '-','Color','#702963','LineWidth',lw);
xlim([-2.7,9]);
title('$\sigma^2$', 'Interpreter', 'latex')
set(gca, 'XTick', [], 'XTickLabel', [],'YTick', [], 'YTickLabel', []);
hold off;
