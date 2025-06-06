function kendallTau = calculate_kendall_tau(yInts, geo)
    % Filter valid rays using the helper function
    validIndices = is_valid_ray(yInts, geo);
    
    % Extract the first column of the valid rays only
    validYInts = yInts(validIndices, 1);
    
    % Calculate Kendall's Tau to measure the order
    if isempty(validYInts) || isscalar(validYInts)
        kendallTau = 0; % Not enough data to calculate Kendall's Tau
    else
        kendallTau = abs(corr(validYInts, (1:length(validYInts))', 'type', 'Kendall'));
    end
end


function validIndices = is_valid_ray(yInts, geo)
    % Helper function to check if a ray is valid (green rays)
    % A valid ray satisfies the following:
    % 1. Its starting y-coordinate is within bounds.
    % 2. It is not zero at the starting y-coordinate.
    
    validIndices = (yInts(:, 1) <= (geo.arrayW + geo.k)) & ...
                   (yInts(:, 1) >= -(geo.arrayW + geo.k)) & ...
                   (yInts(:, 1) ~= 0);
end