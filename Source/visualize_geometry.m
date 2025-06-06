function visualize_geometry(x, f, geo, xInts, yInts, intensityProfile, scale, savePath)
% VISUALIZE_GEOMETRY Visualizes the tank setup, ray paths, intensity profile,
% and displays device geometries (x) and performance metrics (f) on top.
%
%   Inputs:
%       x                - Geometries of the device (array)
%       f                - Performance metrics (array)
%       geo, xInts, yInts- Geometry and ray data for the tank setup.
%       intensityProfile - Vector of intensity values (zeros will be ignored).
%       scale            - Scaling factor for geometry elements.
%       savePath         - (Optional) File path for saving the figure.
%
%   Example:
%       visualize_geometry(x, f, geo, xInts, yInts, intensityProfile, scale, 'tankFig.png');

% Define theta for drawing circles/ellipses.
theta = 0:0.001:2*pi;

% Create a new figure.
figure;

%% Create Left Axes for the Geometry Plot
% Define axes position: [left, bottom, width, height] (normalized units)
leftAxPos = [0.1, 0, 0.6, 1];  % adjust as needed
ax1 = axes('Position', leftAxPos);
hold(ax1, 'on');
axis(ax1, 'square');

% --- Plot Tank Geometry Elements in ax1 ---
% (Ensure that your helper functions plot into the current axes.)
axes(ax1); % make sure ax1 is current
plot_ellipse_and_circles(geo, scale, theta);

% Plot rear lens based on its type.
switch(geo.rearLensType)
    case 'conic'
        bEll2 = geo.rearLensGeo(1);
        ecc2 = geo.rearLensGeo(2);
        aEll2 = sqrt(bEll2^2/(1 + ecc2^2));
        hFace2 = geo.d2 - aEll2;
        kFace2 = 0;
        plot(hFace2 + aEll2 * cos(theta([4713:end, 1:1571])), ...
            kFace2 + bEll2 * sin(theta([4713:end, 1:1571])), ...
            'LineWidth', 1.5, 'Color', 'blue');
    case 'poly'
        yRange = linspace(-100, 100, 1000);
        xPoly = polyval(geo.rearLensGeo, yRange);
        plot(xPoly, yRange, 'LineWidth', 1.5, 'Color', 'blue');
    case 'highPoly'
        yRange = linspace(-100, 100, 1000);
        xPoly = polyval(geo.rearLensGeo, yRange);
        plot(xPoly, yRange, 'LineWidth', 1.5, 'Color', 'blue');
    case 'aspheric'
        yRange = linspace(-100, 100, 1000);
        r = geo.rearLensGeo(1);
        k = geo.rearLensGeo(2);
        A4 = geo.rearLensGeo(3);
        A6 = geo.rearLensGeo(4);
        d2 = geo.rearLensGeo(5);
        asp = @(y) -y^2/(r*(1+sqrt(1 - (1 + k)*y^2/r^2))) - A4*y^4 - A6*y^6 + d2;
        xAsp = zeros(size(yRange));
        for i = 1:length(yRange)
            xAsp(i) = asp(yRange(i));
        end
        plot(xAsp, yRange, 'LineWidth', 1.5, 'Color', 'blue');
end

% Plot detector boundaries and central scatter points.
plot_detectors_and_central_points(geo, scale, xInts);

% Draw ray paths.
plot_ray_paths(xInts, yInts, geo, scale);

% Finalize the geometry plot.
xlim(ax1, [min([-x(1), -100]), max([100, geo.d2+x(6)])]);
ylim(ax1, [-90, 90]);
daspect(ax1, [1, 1, 1]);
hold(ax1, 'off');

%% Create Right Axes for the Intensity Profile Plot
% Define the axes so that they share the same vertical extent.
rightAxPos = [0.8, 0.3, 0.1, 0.4];  % adjust as needed
ax2 = axes('Position', rightAxPos);
title('Intensity Plot', 'Interpreter','latex');
hold(ax2, 'on');

% --- Process the Intensity Profile ---
% Remove zero entries.
nonZeroMask = intensityProfile ~= 0;
idxNonZero = find(nonZeroMask);
intensityNonZero = intensityProfile(nonZeroMask);

% If there is data, fit a polynomial curve.
if ~isempty(idxNonZero)
    % Choose a polynomial degree (up to degree 5 or fewer if not enough points).
    %degree = min(16, length(idxNonZero)-1);
    %p = polyfit(idxNonZero, intensityNonZero, degree);
    % Evaluate the fit on a smooth grid.
    %idxFit = linspace(min(idxNonZero), max(idxNonZero), 200);
    %intensityFit = polyval(p, idxFit);

    % --- Plot the Data ---
    % Note: To "rotate" the plot so that the index is vertical, we swap x and y.
    % Plot raw data as small black circles.
    plot(ax2, intensityNonZero, idxNonZero, 'ko', 'MarkerSize', 3, 'DisplayName', 'Raw Data');
    % Optionally, plot the fitted curve (uncomment if desired).
    %plot(ax2, intensityFit, idxFit, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Fitted Curve');
end
xlim(ax2, [0, 200]);
ylim(ax2, [0, 2048]);
hold(ax2, 'off');

n = numel(x);

% 1) Build the LaTeX pieces for headers and values
headersLatex = arrayfun(@(i) sprintf('x_{%d}', i), 1:n, 'UniformOutput', false);
valuesLatex  = arrayfun(@(v) sprintf('%.4f',   v), x,   'UniformOutput', false);

% 2) Make a column specifier of “n” centered columns
colSpec = repmat('c', 1, n);  

% 3) Join into a single array environment
arrayStr = [ ...
  '$\displaystyle\begin{array}{' colSpec '}'   ...  % start math mode & array
     strjoin(headersLatex, ' & ') '\\'          ...  % first row: headers
     strjoin(valuesLatex,  ' & ')                ...  % second row: values
  '\end{array}$'                                  ...  % end array & math mode
];

% 4) Third line as its own LaTeX math string
rhoStr = sprintf('\n$\\rho_{eff}=%.2f\\;,\\;|\\tau|=%.2f\\;,\\;CV_{I}=%.2f$', ...
                 f(1), f(2), f(3));

% 5) Draw the textbox with the LaTeX interpreter
annotation('textbox', [0.1, 0.86, 0.8, 0.12], ...
    'String',             { arrayStr; rhoStr }, ...
    'Interpreter',        'latex', ...
    'FontSize',           10, ...
    'HorizontalAlignment','center', ...
    'LineStyle',          'none');



%% Save the Figure (if a save path is provided)
if nargin >= 8 && ~isempty(savePath)
    saveas(gcf, savePath);
    disp(['Figure saved to: ', savePath]);
end

% Optionally close the figure to prevent clutter.
close(gcf);
end
