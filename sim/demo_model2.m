function demo_model2(posterior, center_axis, estx, conf_radius, conf_radius_axis, costFun, erCDF, gainFun)

lw = 1.5;

%% final posterior
f1 = figure;
set(gcf, 'Position',[0 0 500 200])
set(gca, 'LineWidth', 1, 'FontSize', 15,'TickDir', 'out'); hold on
plot(center_axis, posterior, 'k','LineWidth', lw);
l1 = xline(estx, 'r', 'LineWidth', lw);

line1 = estx - conf_radius;
line2 = estx + conf_radius;

% Get the y-values (posterior values) at line1 and line2 positions
y_line1 = interp1(center_axis, posterior, line1);
y_line2 = interp1(center_axis, posterior, line2);

% Plot the two symmetric lines up to the intersection with the posterior curve
lc = plot([line1, line1], [0, y_line1], 'b-', 'LineWidth', 1.5);  % Left line
uc = plot([line2, line2], [0, y_line2], 'b-', 'LineWidth', 1.5);  % Right line

% Find the region between the two lines for shading
x_shaded = center_axis(center_axis >= line1 & center_axis <= line2);
y_shaded = posterior(center_axis >= line1 & center_axis <= line2);

% Shade the region between the two lines under the curve
fill([x_shaded, fliplr(x_shaded)], [y_shaded, zeros(size(y_shaded))], 'k', 'FaceAlpha', 0.3);

% Add labels
xlim([-60, 60])
xlabel('Stimulus location (cm)');
ylabel('Posterior probability');
legend([l1, lc, uc], {'Final estimate','Lower confidence bound','Higher confidence bound'})
legend boxoff
set(f1, 'DefaultFigureRenderer', 'painters');
saveas(f1, 'demo1-1', 'epsc')

%% probability mass within the confidence radius

f2 = figure;
set(gcf, 'Position',[0 0 500 400])

subplot(3,1,1)
set(gca, 'LineWidth', 1, 'FontSize', 15,'TickDir', 'out'); hold on
plot(conf_radius_axis, erCDF, 'k','LineWidth',lw);
ylabel('Probability mass')
xline(conf_radius,'b-','LineWidth',lw);
xlim([0, 50])

%% cost function
subplot(3,1,2)
set(gca, 'LineWidth', 1, 'FontSize', 15,'TickDir', 'out'); hold on
plot(conf_radius_axis, costFun,'k', 'LineWidth',lw);
ylabel('Points')
xline(conf_radius,'b-','LineWidth',lw);
xlim([0, 50])

%% expected gain
subplot(3,1,3);
set(gca, 'LineWidth', 1, 'FontSize', 15,'TickDir', 'out'); hold on
plot(conf_radius_axis, gainFun, 'k','LineWidth',lw);
xlabel('Confidence radius (cm)');
ylabel('Expected gain')
xline(conf_radius,'b-','LineWidth',lw);
xlim([0, 50])
ylim([0,1])

set(f2, 'DefaultFigureRenderer', 'painters');
saveas(f2, 'demo1-2', 'epsc')

end