%% Define Parameters for Data and Plotting
% Adjust these as needed
gelIndex = 1.4225; %Flexidose3D
% gelIndex = 1.34468;%ClearView
% gelIndex = 1.3319; %Water
rearLensType = 'highPoly';

% Define the data directory based on the rearLensType.
% For example, if rearLensType = 'aspheric', the CSV is expected at:
%   Data/aspheric/ParetoResults.csv
dataDir = fullfile(pwd, 'Data', 'HighPoly');
csvFile = fullfile(dataDir, 'ParetoResults.csv');

%% Read ParetoResults.csv and Filter Solutions with rayCount = 10000
% Read the CSV file from the specified data directory
resultsTable = readtable(csvFile);

% Filter rows for which the rayCount equals 10000
rows = resultsTable(resultsTable.rayCount == 100000, :);
rows = sortrows(rows, 'Objective1');

% Check if any solution with 10000 rays is available
if isempty(rows)
    error('No solution with rayCount = 100000 found in %s.', csvFile);
end

% Determine the number of solutions to process
numSolutions = height(rows);
fprintf('Found %d solution(s) with rayCount = 100000.\n', numSolutions);

%% Prepare the Output Directory for Images
% Create an output directory "img" under the data folder if it does not exist
imgDir = fullfile(dataDir, 'img');
if ~exist(imgDir, 'dir')
    mkdir(imgDir);
    fprintf('Created image directory: %s\n', imgDir);
end

%% Loop Over Each Solution and Plot the Geometry
for i = 1:numSolutions
    % Extract the i-th row from the filtered results
    selectedRow = rows(i, :);
    
    % Extract the design vector x (columns x1 to x10)
    x = table2array(selectedRow(:, 1:14));
    
    % Define a unique save path for each plot inside the "img" folder
    savePath = fullfile(imgDir, sprintf('SolidTankGeoPlot_%d.png', i));
    
    % Call the SolidTankPlot function to plot the geometry
    % SolidTankPlot should be defined to accept parameters (x, gelIndex, rearLensType, savePath)
    SolidTankPlot(x, gelIndex, rearLensType, savePath);
    
    % Display a message indicating that the plot was saved
    fprintf('Plot %d saved to: %s\n', i, savePath);
end
