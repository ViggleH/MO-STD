function stdIntensity = calculate_cv_intensity(intensityProfile)
if mean(intensityProfile) ~= 0
    stdIntensity = std(intensityProfile)/mean(intensityProfile); 
else
    stdIntensity = 10e6;
end
