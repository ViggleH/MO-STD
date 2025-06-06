% List of rear‐lens types
rearLensTypes = {'Conic', 'Conic-Single', 'Poly', 'Poly-single', 'Aspheric', 'Aspheric-Single'};

% Preallocate cell arrays
ParetoFronts      = cell(size(rearLensTypes));
SelectedSolutions = cell(size(rearLensTypes));
Hypervolumes      = zeros(size(rearLensTypes));

%% 1) Read & filter Pareto fronts (rayCount == 100000), then sort by Objective1
for i = 1:numel(rearLensTypes)
    % Construct path to the CSV
    csvPath = fullfile('.', 'Data', rearLensTypes{i}, 'ParetoResults.csv');
    
    % Read the CSV into a table
    resultsTable = readtable(csvPath);
    
    % Filter rows where rayCount == 100000
    rows100k = resultsTable(resultsTable.rayCount == 100000, :);
    
    % Sort those rows by 'Objective1'
    ParetoFronts{i} = sortrows(rows100k, 'Objective1');
end

%% 2) Select only those solutions with Objective3 <= -0.99
for i = 1:numel(rearLensTypes)
    PF = ParetoFronts{i};
    if ~isempty(PF)
        SelectedSolutions{i} = PF(PF.Objective3 <= -0.99, :);
    else
        SelectedSolutions{i} = PF;  % empty table
    end
end

%% 3) Compute the global Ideal (min) and Nadir (max) across all selectedSolutions
ideal1 =  Inf; ideal2 =  Inf; ideal3 =  Inf;
nadir1 = -Inf; nadir2 = -Inf; nadir3 = -Inf;

for i = 1:numel(rearLensTypes)
    T = SelectedSolutions{i};
    if ~isempty(T)
        ideal1 = min(ideal1, min(T.Objective1));
        ideal2 = min(ideal2, min(T.Objective2));
        ideal3 = min(ideal3, min(T.Objective3));
        
        nadir1 = max(nadir1, max(T.Objective1));
        nadir2 = max(nadir2, max(T.Objective2));
        nadir3 = max(nadir3, max(T.Objective3));
    end
end

% Check if any solutions were found at all
if isinf(ideal1)
    error('No selected solutions found across any rear-lens type.');
end

Ideal = [ideal1, ideal2, ideal3];
Nadir = [nadir1, nadir2, nadir3];

fprintf('Global Ideal  = [%.4f, %.4f, %.4f]\n', Ideal);
fprintf('Global Nadir  = [%.4f, %.4f, %.4f]\n', Nadir);

%% 4) Compute 3D hypervolume for each rear‐lens type
for i = 1:numel(rearLensTypes)
    T = SelectedSolutions{i};
    if isempty(T)
        Hypervolumes(i) = 0;
    else
        % Extract the 3‐column matrix of objectives [f1, f2, f3]
        PF_mat = T{:, {'Objective1','Objective2','Objective3'}};
        Hypervolumes(i) = computeHypervolume3D(PF_mat, Nadir, Ideal);
    end
    fprintf('Hypervolume (%s) = %.6f\n', rearLensTypes{i}, Hypervolumes(i));
end