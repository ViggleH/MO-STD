function plot_ellipse_and_circles(geo, scale, theta)
    % Plot front ellipse
    hFace = -geo.d1 * scale + (geo.bEll^2 / (1 + geo.ecc^2))^0.5;
    kFace = 0;
    plot(hFace + (geo.bEll^2 / (1 + geo.ecc^2))^0.5 * cos(theta(1571:4713)), ...
         kFace + geo.bEll * sin(theta(1571:4713)), 'LineWidth', 1.5, 'Color', 'blue');

 
    % Plot circular boundaries (gel, container, bore radii)
    plot(geo.k + geo.r1 * sin(theta), geo.r1 * cos(theta), 'LineWidth', 1, 'Color', 'blue');
    plot(geo.k + geo.r2 * sin(theta), geo.r2 * cos(theta), 'LineWidth', 1, 'Color', 'magenta');
    plot(geo.k + geo.r3 * sin(theta), geo.r3 * cos(theta), 'LineWidth', 1, 'Color', 'black');
end