% Parameters for the narrow Gaussian
mu1 = 0; % Mean
sigma1 = 0.5; % Standard deviation

% Parameters for the wider Gaussian
mu2 = 3; % Mean
sigma2 = 1.5; % Standard deviation

% Weights for the two Gaussians
w1 = 0.4; % Weight for the narrow Gaussian
w2 = 0.6; % Weight for the wider Gaussian

% Define the x-axis range
x = linspace(-5, 7, 1000);

% Compute the two Gaussian distributions
gauss1 = w1 * normpdf(x, mu1, sigma1);
gauss2 = w2 * normpdf(x, mu2, sigma2);

% Sum of the weighted Gaussians
combined = gauss1 + gauss2;
norm_combined = combined./sum(combined);

% Plot the combined Gaussian distribution
figure;
plot(x, norm_combined, 'k', 'LineWidth', 2);
hold on;

% Customize the plot
xlabel('x');
ylabel('Probability Density');
title('Weighted Combination of Two Gaussian Distributions');
grid on;
hold off;