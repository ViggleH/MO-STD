function plot_detectors_and_central_points(geo, scale, xInts)
    % Plot detector array
    line([geo.det * scale, geo.det * scale], ...
         [(geo.arrayW + geo.k) * scale, -(geo.arrayW + geo.k) * scale], ...
         'LineWidth', 5, 'Color', 'black');

    % Plot leftmost intersection line for context
    line([min(xInts(:)) * scale, min(xInts(:)) * scale], ...
         [5 * scale, -5 * scale], 'LineWidth', 1, 'Color', 'black');
    line([min(xInts(:)) * scale, (min(xInts(:)) - 10) * scale], ...
         [5 * scale, 5 * scale], 'LineWidth', 2, 'Color', 'black');
    line([min(xInts(:)) * scale, (min(xInts(:)) - 10) * scale], ...
         [-5 * scale, -5 * scale], 'LineWidth', 2, 'Color', 'black');

    % Scatter key central points
    scatter(0, 0, 'black', 'filled'); % Tank center
    scatter(geo.k, 0, 'red', 'filled'); % Bore center
end