% Multifidelity Bi-objective Optimization using Genetic Algorithm
% Maximizing Both Objectives

% Clear workspace and command window
clear; clc; close all;
diary('.\Data\Aspheric-Single\output.txt');
t_total = tic;
%% Define Problem Parameters
gelIndex = 1.4225; %Flexidose3D
% gelIndex = 1.34468;%ClearView
% gelIndex = 1.3319; %Water
rearLensType = 'aspheric';

nVars = 10; % Number of decision variables (adjust as needed)
lb = [200, -40, 40, 40, 0, 10, 0, -10, -10, -10]; % Lower bounds for decision variables
ub = [500, 40, 100, 160, 5, 400, 500, 10, 10, 10]; % Upper bounds for decision variables

% Define rayCount levels (from low to high fidelity)
rayCounts = [100000]; % Example rayCount levels

% Number of Pareto solutions to store
popSizes = [256]*8;
% Maximum number of generations for the genetic algorithm
generations = [128];

% Initialize cell arrays to store Pareto fronts
paretoFronts = cell(length(rayCounts),1);
paretoX = cell(length(rayCounts),1);

% Initialize a table to store all results
% Preallocate with variable names
resultTable = table();

%% Optimization Options
% Common options for gamultiobj
options = optimoptions('gamultiobj');
options.Display = 'iter'; % Display output at each generation
%options.PlotFcn = {@gaplotpareto}; % Optional: Plot Pareto front at each iteration
options.UseParallel = true;
options.UseVectorized = false;
options.MaxTime = 2994.80;


% Optional: You can further customize the genetic algorithm options as needed
% For example:
% options.CrossoverFraction = 0.8;
% options.MutationRate = 0.1;
% options.EliteCount = round(0.05 * popSize);
% options.SelectionFcn = @selectiontournament;

%% Iterate Over RayCount Levels
for i = 1:length(rayCounts)
    options.MaxGenerations = generations(i);
    options.PopulationSize = popSizes(i); % Size of the population
    currentRayCount = rayCounts(i);
    fprintf('Starting optimization with rayCount = %d\n', currentRayCount);
    
    % Define the objective function handle with current rayCount
    % Since gamultiobj minimizes, we negate the objectives to maximize
    objFun = @(x) TestObj(x, currentRayCount, gelIndex, rearLensType);
    
    % Set Initial Population
    if i == 1
        % For the first level, use a random initial population
        initialPop = repmat(lb, popSizes(i), 1) + rand(popSizes(i), nVars) .* repmat((ub - lb), popSizes(i), 1);
    else
        % For subsequent levels, use the previous Pareto-optimal solutions with best Objective1 performance.
        % Retrieve previous decision vectors and objective values.
        prevX = paretoX{i-1};
        prevF = paretoFronts{i-1};
        
        % Sort the previous solutions based on Objective1 in ascending order.
        % (Assuming lower Objective1 is better; if not, change 'ascend' to 'descend')
        [~, sortIdx] = sort(prevF(:,1), 'ascend');
        prevX_sorted = prevX(sortIdx, :);
        numPrev = size(prevX_sorted, 1);
        
        if numPrev >= popSizes(i)
            % Use the best popSizes(i) solutions
            initialPop = prevX_sorted(1:popSizes(i), :);
        else
            % If fewer, use all and fill the remainder with random points.
            numRandom = popSizes(i) - numPrev;
            randomPop = repmat(lb, numRandom, 1) + rand(numRandom, nVars) .* repmat((ub - lb), numRandom, 1);
            initialPop = [prevX_sorted; randomPop];
        end
    end
    
    % Update options with the chosen Initial Population
    options.InitialPopulationMatrix = initialPop;
    
    % Run the multi-objective genetic algorithm
    tic;
    [xPareto, fPareto] = gamultiobj(objFun, nVars, [], [], [], [], lb, ub, options);
    elapsedTime = toc;
    
    % Store the Pareto fronts and solutions
    paretoX{i} = xPareto;
    paretoFronts{i} = fPareto;
    
    fprintf('Optimization with rayCount = %d completed. Number of Pareto solutions: %d\n\n', currentRayCount, size(fPareto,1));
    
    % Create a temporary table for current rayCount
    tempTable = array2table([xPareto, repmat(currentRayCount, size(xPareto,1),1), fPareto], ...
        'VariableNames', {'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9', 'x10', 'rayCount', 'Objective1', 'Objective2', 'Objective3'});
    fprintf('Computation Time: %.2f seconds\n\n', elapsedTime);
    
    % Append additional data (e.g., computing time) to the temporary table
    tempTable = array2table([xPareto, ...
        repmat(currentRayCount, size(xPareto,1),1), ...
        fPareto, ...
        repmat(elapsedTime, size(xPareto,1),1)], ...
        'VariableNames', {'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9', 'x10', 'rayCount', 'Objective1', 'Objective2', 'Objective3', 'ComputingTime'});
    
    % Append to the master results table
    resultTable = [resultTable; tempTable];
end


%% Plot All Pareto Fronts and Save Figures
% Create a new figure for the Pareto front plot
hFig = figure('Units','pixels',...
       'Position',[100 100 640 160]);
hold on;
colors = lines(length(rayCounts));  % Distinct colors for each fidelity level

for i = 1:length(rayCounts)
    % Find the indices for points meeting the condition for group i
    validIdx = paretoFronts{i}(:,3) <= -0.99;
    
    % Only call scatter if there are any valid indices
    if any(validIdx)
        scatter(paretoFronts{i}(validIdx, 1), paretoFronts{i}(validIdx, 2), ...
            36, colors(i,:), 'DisplayName', sprintf('$N_{\\mathrm{rays}}$ = %d', rayCounts(i)));
    end
end

xlabel('$-\rho_{\mathrm{eff}}$', 'Interpreter', 'latex');
ylabel('$\mathrm{CV_{I}}$', 'Interpreter', 'latex');
axScaled = findobj(hFig, 'Type', 'axes');
set(axScaled, 'XLim', [-1, -0.8], 'YLim', [0, 20]);

title('Pareto Fronts at Different $N_{\mathrm{rays}}$ (Aspheric Lens, FlexyDos3D)', Interpreter='latex');
legend('show', 'Interpreter', 'latex', 'Location', 'best');
grid on;
hold off;
%% Define the Output Folder
% Set the output folder (for example, 'Data/Conic')
outputFolder = fullfile('Data', 'Aspheric-Single');

% Create the folder if it does not exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Save the Figure as a .fig File
figFile = fullfile(outputFolder, 'ParetoFronts.fig');
savefig(hFig, figFile);
fprintf('Figure saved as .fig: %s\n', figFile);

%% Save a Scaled PNG Version of the Figure
% To preserve the original figure, create a copy of it.
hFigScaled = copyobj(hFig, 0);  % Copy to the root (0) to create a new independent figure

% Find the axes in the copied figure and adjust axis limits
axScaled = findobj(hFigScaled, 'Type', 'axes');
set(axScaled, 'XLim', [-2, -1.8], 'YLim', [0, 2]);

% Save the scaled figure as a PNG file
pngFile = fullfile(outputFolder, 'ParetoFronts_scaled.png');
saveas(hFigScaled, pngFile);
fprintf('Scaled PNG saved: %s\n', pngFile);

% Close the scaled figure copy to free memory
close(hFigScaled);

%% Save the Master Table and Additional Results
% Save the master table to a .mat file
matFile = fullfile(outputFolder, 'ParetoResults.mat');
save(matFile, 'resultTable', 'paretoX', 'paretoFronts');
fprintf('Master MAT file saved: %s\n', matFile);

% Save the table as a CSV file
csvFile = fullfile(outputFolder, 'ParetoResults.csv');
writetable(resultTable, csvFile);
fprintf('Master CSV file saved: %s\n', csvFile);

fprintf('All results have been successfully saved.\n');
toc(t_total)
diary off;
%% Example Objective Function
% Replace this with your actual SolidTankObj function
function F = TestObj(X, rayCount, gelIndex, rearLensType)
% TestObj Vectorized multi-objective wrapper around SolidTankObj.
%   X            — P×nVars matrix, each row is one candidate design
%   rayCount     — scalar fidelity parameter
%   gelIndex     — scalar or index parameter for the gel
%   rearLensType — scalar or code for the rear-lens type
%
%   F            — P×3 matrix of objective values, where
%                  F(:,1) = -effRad, F(:,2) = CV, F(:,3) = -tau

% Number of candidate solutions
P = size(X, 1);
% Preallocate output
F = zeros(P, 3);

% Loop over each candidate
for k = 1:P
    xk = X(k, :);
    try
        % Evaluate the true objectives
        [effRad, CV, tau] = SolidTankObj(xk, rayCount, gelIndex, rearLensType);

        % Negate where needed
        obj = [ -effRad, CV, -tau ];

        % Check for non‑real / invalid numbers
        if any(~isreal(obj)) || any(isnan(obj)) || any(isinf(obj))
            error('TestObj:InvalidOutput', ...
                'Non‑real or invalid objective for candidate %s', mat2str(xk));
        end

        % Store the valid objectives
        F(k, :) = obj;

    catch ME
        % Properly pass identifier + formatted message:
        warning(ME.identifier, ...
            'TestObj failure: %s\nInput x: %s\nReturning penalty values.', ...
            ME.message, mat2str(xk));
        F(k,:) = [1e6,1e6,1e6];
    end

end
end
