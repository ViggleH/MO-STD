function score3 = calculate_effective_radius(geo, xInts, yInts)
% Identify the shortest ray that passes through the center of the tank
shortRay = zeros(size(xInts, 1), 1);
for i = 1:size(xInts, 1)
    if is_valid_ray(yInts(i, :), geo)
        % Compute chord length
        intPts = [xInts(i, 6), yInts(i, 6); xInts(i, 5), yInts(i, 5)];
        shortRay(i) = pdist(intPts, 'euclidean');
    else
        shortRay(i) = NaN;
    end
end

% Find the effective radius
if isempty(shortRay) || all(isnan(shortRay))
    effRad = 0; % No valid rays
else
    % Identify continuous groups of non-NaN values
    isValid = ~isnan(shortRay); % Logical array indicating non-NaN values
    groups = bwlabel(isValid); % Label continuous groups of non-NaN values

    % Find the largest continuous group
    groupSizes = histcounts(groups, 'BinLimits', [1, max(groups)], 'BinWidth', 1);
    [~, largestGroup] = max(groupSizes);

    % Mask to retain only the largest group
    largestGroupMask = (groups == largestGroup);

    % Filter valid rays from the largest group
    validRays = shortRay(largestGroupMask);
    nonNanIndices = find(largestGroupMask);

    % If there are no valid rays after filtering
    if isempty(validRays)
        effRad = 0; % No valid rays in the largest continuous group
    else
        % Find the index of the ray closest to the center (counted top-down)
        [~, furthestIndex] = min(abs(validRays - validRays(1))); % Closest to center
        effRadSelect = nonNanIndices(furthestIndex);

        % Calculate the chord for the furthest ray
        chordx = [xInts(effRadSelect, 6), xInts(effRadSelect, 5)];
        chordy = [yInts(effRadSelect, 6), yInts(effRadSelect, 5)];

        % Math for finding crossing point between center and line segment
        slope = diff(chordy) / diff(chordx); % Correct slope calculation
        yIntercept = chordy(1) - slope * chordx(1); % y-intercept of chord
        yIntercept2 = geo.k - (-1 / slope) * geo.realh; % y-intercept of line from center
        xIntersect = (yIntercept2 - yIntercept) / (slope - (-1 / slope));
        yIntersect = slope * xIntersect + yIntercept;

        % Calculate effective radius (distance from center to intersection)
        effRad = norm([xIntersect - geo.realh, yIntersect - geo.k]);
    end
end


% Compute the score
%score3 = 0.5 * tanh((10) * (effRad / geo.r1 - 0.96) * 2 * pi) + 0.5;
score3 = effRad / geo.r1;

end

function valid = is_valid_ray(yInts, geo)
    % Helper function to check if a ray is valid (green rays)
    % A valid ray satisfies the following:
    % 1. Its starting y-coordinate is within bounds.
    % 2. It is not zero at the starting y-coordinate.
    
    valid = all(yInts(:, 1) <= (geo.arrayW + geo.k)) && ...
            all(yInts(:, 1) >= -(geo.arrayW + geo.k)) && ...
            all(yInts(:, 1) ~= 0);
end

%{
function valid = is_valid_ray(xInts, yInts, geo)
% Helper function to check if a ray is valid
valid = (max(yInts) < (geo.arrayW + geo.k)) && ...
    (min(yInts) > -(geo.arrayW + geo.k)) && ...
    (nnz(xInts) > 9) && ...
    (xInts(1) == geo.det);
end
%}