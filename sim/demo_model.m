function demo_model(aloc, vloc, tt, fixP, post, shat, conf_radius, conf_radius_axis, costFun, erCDF, gainFun)

lw = 1.5;

% Figure handle (set your desired figure number)
figHandle = findobj('Type', 'figure', 'Name', 'demo');

if isempty(figHandle)
    % If the figure does not exist, create a new one
    figure('Name', 'demo');
    set(gcf, 'Position',[0 00 500 800])
    hold on;

    %% final posterior
    subplot(4,1,1);
    set(gca, 'LineWidth', 1, 'FontSize', 15,'TickDir', 'out'); hold on
    posterior = squeeze(post(aloc, vloc, 1, 1, tt, :));
    plot(fixP.center_axis, posterior, 'k','LineWidth', lw);
    shat_value = shat(aloc, vloc, 1, 1, tt);
    l1 = xline(shat_value, 'r', 'LineWidth', lw);

    line1 = shat_value - conf_radius;
    line2 = shat_value + conf_radius;

    % Get the y-values (posterior values) at line1 and line2 positions
    y_line1 = interp1(fixP.center_axis, squeeze(post(aloc, vloc, 1, 1, tt, :)), line1);
    y_line2 = interp1(fixP.center_axis, squeeze(post(aloc, vloc, 1, 1, tt, :)), line2);

    % Plot the two symmetric lines up to the intersection with the posterior curve
    lc = plot([line1, line1], [0, y_line1], 'b-', 'LineWidth', 1.5);  % Left line
    uc = plot([line2, line2], [0, y_line2], 'b-', 'LineWidth', 1.5);  % Right line

    % Find the region between the two lines for shading
    x_shaded = fixP.center_axis(fixP.center_axis >= line1 & fixP.center_axis <= line2);
    y_shaded = squeeze(post(aloc, vloc, 1, 1, tt, fixP.center_axis >= line1 & fixP.center_axis <= line2));

    % Shade the region between the two lines under the curve
    fill([x_shaded, fliplr(x_shaded)], [y_shaded', zeros(size(y_shaded'))], 'k', 'FaceAlpha', 0.3);

    % Add labels
    xlim([-60, 60])
    xlabel('Stimulus location (cm)');
    ylabel('Posterior probability');
    legend([l1, lc, uc], {'Final estimate','Lower confidence bound','Higher confidence bound'})
    legend boxoff

    hold off;

else

    % If the figure exists, make it the current figure
    figure(figHandle);
    hold on;
    disp('Adding to the existing figure');

    %% probability mass within the confidence radius
    subplot(4,1,2);
    set(gca, 'LineWidth', 1, 'FontSize', 15,'TickDir', 'out'); hold on
    plot(conf_radius_axis, erCDF, 'k','LineWidth',lw);
    xlabel('Confidence radius (cm)');
    ylabel('Probability mass')
    xline(conf_radius,'b-','LineWidth',lw);
    xlim([0, 50])

    %% cost function
    subplot(4,1,3);
    set(gca, 'LineWidth', 1, 'FontSize', 15,'TickDir', 'out'); hold on
    plot(conf_radius_axis, costFun,'k', 'LineWidth',lw);
    xlabel('Confidence radius (cm)');
    ylabel('Points')
    xline(conf_radius,'b-','LineWidth',lw);
    xlim([0, 50])

    %% expected gain
    subplot(4,1,4);
    set(gca, 'LineWidth', 1, 'FontSize', 15,'TickDir', 'out'); hold on
    plot(conf_radius_axis, gainFun, 'k','LineWidth',lw);
    xlabel('Confidence radius (cm)');
    ylabel('Expected gain')
    xline(conf_radius,'b-','LineWidth',lw);
    xlim([0, 50])
    ylim([0,1])

    set(gcf, 'DefaultFigureRenderer', 'painters');
    saveas(gcf, 'demo2', 'epsc')

end

end