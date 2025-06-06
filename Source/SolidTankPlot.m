function SolidTankPlot(x, gelIndex, rearLensType, savePath)
    rayCount = 50;
    scale = 1;
    laserType = 'fan';
    N = 512;
    numDet = 2048; % Number of detectors
    
    % Step 1: Setup geometry and refractive parameters
    geo = setup_geometry(rayCount, numDet, laserType, N, x, gelIndex, rearLensType);
    geo2 = setup_geometry(100000, numDet, laserType, N, x, gelIndex, rearLensType);
    
    % Step 2: Compute intersection points and intensities
    [xInts, yInts, intensityProfile] = compute_intersections(geo);
    [xInts2, yInts2, intensityProfile2] = compute_intersections(geo2);

    effRad = calculate_effective_radius(geo2, xInts2, yInts2);
    tau =  calculate_kendall_tau(yInts2,geo2);
    intensityCV = calculate_cv_intensity(intensityProfile2); 
    
    % Step 3: Visualize the setup and ray paths
    visualize_geometry(x, [effRad, tau, intensityCV], geo, xInts, yInts, intensityProfile2, scale, savePath)
    
end