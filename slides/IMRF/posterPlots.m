% Define parameters for the first Gaussian distribution
mu1 = 0; % mean
sigma1 = 1; % standard deviation

% Define parameters for the second Gaussian distribution
mu2 = 5; % mean
sigma2 = 1.5; % standard deviation

% Create the x range for the PDF
x = linspace(mu1 - 4*sigma1, mu2 + 4*sigma2, 1000);

% Calculate the PDFs
pdf1 = normpdf(x, mu1, sigma1);
pdf2 = normpdf(x, mu2, sigma2);

% Calculate reliability (inverse of variance)
reliability1 = 1 / sigma1^2;
reliability2 = 1 / sigma2^2;

% Weighted combination of the two PDFs
combined_mean = (reliability1 * mu1 + reliability2 * mu2) / (reliability1 + reliability2);
combined_variance = 1 / (reliability1 + reliability2);
combined_pdf = normpdf(x, combined_mean, sqrt(combined_variance));

% Plot the PDFs
figure;

subplot(3,1,1);
% Plot the grey shade under the first Gaussian curve
fill([x, fliplr(x)], [combined_pdf, zeros(size(combined_pdf))], [0.8 0.8 0.8], 'EdgeColor', 'none');
hold on;
lw = 2;
% Plot the Gaussian curves
plot(x, pdf1, '-','Color','#B4222D','LineWidth',lw);
plot(x, pdf2, '-','Color','#2270B0','LineWidth',lw);
plot(x, combined_pdf, 'k-','LineWidth',lw);

% Add labels and legend
legend('','${m_V}$','${m_A}$','$\hat{s}$', 'Interpreter', 'latex');
xlabel('Value');
ylabel('Probability Density');
title('Gaussian Distributions and Combined Distribution');
hold off;

subplot(3,1,2);
% Plot the grey shade under the first Gaussian curve
fill([x, fliplr(x)], [pdf1, zeros(size(pdf1))], [0.8 0.8 0.8], 'EdgeColor', 'none');
hold on;
lw = 2;
% Plot the Gaussian curves
plot(x, pdf1, '-','Color','#B4222D','LineWidth',lw);
plot(x, pdf2, '-','Color','#2270B0','LineWidth',lw);
plot(x, combined_pdf, 'k-','LineWidth',lw);

% Add labels and legend
legend('','${m_V}$','${m_A}$','$\hat{s}$', 'Interpreter', 'latex');
xlabel('Value');
ylabel('Probability Density');
title('Gaussian Distributions and Combined Distribution');
hold off;
