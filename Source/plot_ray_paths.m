function plot_ray_paths(xInts, yInts, geo, scale)
    for i = 1:size(xInts, 1)
        for j = 1:size(xInts, 2) - 1
            % Define ray endpoints
            xStart = xInts(i, j);
            xEnd = xInts(i, j + 1);
            yStart = yInts(i, j);
            yEnd = yInts(i, j + 1);

            % Determine ray color
            if yInts(i, 1) > (geo.arrayW + geo.k) * scale || ...
               yInts(i, 1) < -(geo.arrayW + geo.k) * scale || ...
               yInts(i, 1) == 0
                rayColor = 'red'; % Invalid ray
            else
                rayColor = [0, 0.5, 0]; % Valid ray (green)
            end

            % Plot ray
            if yStart ~= 0
            line([xStart, xEnd], [yStart, yEnd], 'Color', rayColor);
            end
        end
    end
end