function [effRad, intensityCV, tau] = SolidTankObj(x, rayCount, gelIndex, rearLensType)
    laserType = 'fan';
    N = 512;
    numDet = 2048; % Number of detectors
    
    % Step 1: Setup geometry and refractive parameters
    geo = setup_geometry(rayCount, numDet, laserType, N, x, gelIndex, rearLensType);
    
    % Step 2: Compute intersection points and intensities
    [xInts, yInts, intensityProfile] = compute_intersections(geo);
    
    % Step 3: Visualize the setup and ray paths
    %visualize_geometry(geo, xInts, yInts, scale, savePath);
    
    % Step 4: Calculate scores
    effRad = calculate_effective_radius(geo, xInts, yInts) ; %Effective Radius
    intensityCV = calculate_cv_intensity(intensityProfile); %Inverse Coefficient of Variation of the Intensity Profile
    tau = calculate_kendall_tau(yInts,geo);
    
end