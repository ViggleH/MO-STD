function [xint, yint, aInc, discrim] = LineAspIntersect(xIntercept, yIntercept, slope, aspCoeffs)
% LineAspIntersect computes the intersection of a ray with an aspherical lens.
%
%   [xint, yint, aInc, discrim] = LineAspIntersect(xIntercept, yIntercept, slope, aspCoeffs)
%
% Inputs:
%   xIntercept, yIntercept, slope : define the ray:
%       y = slope*(x - xIntercept) + yIntercept
%   aspCoeffs : a vector containing the aspherical coefficients:
%       aspCoeffs(1) = r   (radius of curvature)
%       aspCoeffs(2) = k   (conic constant)
%       aspCoeffs(3) = A4  (4th-order aspheric coefficient)
%       aspCoeffs(4) = A6  (6th-order aspheric coefficient)
%       aspCoeffs(5) = d2  (offset term)
%
% Outputs:
%   xint, yint : Intersection point between the ray and the aspherical lens.
%   aInc       : Angle of incidence (radians) between the ray and the lens normal.
%   discrim    : 1 if a valid real intersection is found; -1 if no valid intersection.
%
% Improvements:
%   - Uses a custom Newton–Raphson solver.
%   - Checks intermediate and final solutions for non-real values and returns NaN if found.

    % Extract aspherical coefficients.
    r = aspCoeffs(1);
    k = aspCoeffs(2);
    A4 = aspCoeffs(3);
    A6 = aspCoeffs(4);
    d2 = aspCoeffs(5);

    % Define the ray function (explicit in y): x = (y - yIntercept)/slope + xIntercept.
    rayEq = @(y) (y - yIntercept)/slope + xIntercept;

    % Define the function whose root we seek:
    %   F(y) = asph(y) - rayEq(y)
    F = @(y) asph(y, r, k, A4, A6, d2) - rayEq(y);
    % Its derivative:
    dF = @(y) dAsph(y, r, k, A4, A6) - 1/slope;

    % Use Newton–Raphson to solve F(y)=0.
    tol = 1e-8;
    maxIter = 20;
    y_current = yIntercept;  % Use yIntercept as initial guess.
    validSolution = true;
    
    for iter = 1:maxIter
        f_val = F(y_current);
        df_val = dF(y_current);
        
        % If the derivative is zero or nearly zero, break.
        if abs(df_val) < eps
            validSolution = false;
            break;
        end
        
        y_new = y_current - f_val/df_val;
        
        % Check if the new iterate is complex.
        if ~isreal(y_new)
            validSolution = false;
            break;
        end
        
        if abs(y_new - y_current) < tol
            y_current = y_new;
            break;
        end
        
        y_current = y_new;
    end

    % If the solution did not converge or is complex, return invalid.
    if ~validSolution || ~isreal(y_current) || abs(F(y_current)) > tol
        yint = NaN;
        xint = NaN;
        aInc = NaN;
        discrim = -1;
        return;
    else
        yint = y_current;
    end

    % Compute the corresponding x-coordinate on the aspherical surface.
    xint = asph(yint, r, k, A4, A6, d2);
    
    % Simple validity test (application dependent): require xint > 0.
    if xint > 0
        discrim = 1;
    else
        xint = NaN;
        discrim = -1;
    end

    % --- Compute the normal and the incidence angle ---
    % Compute the derivative of the aspherical function at yint.
    dAsp_val = dAsph(yint, r, k, A4, A6);
    
    % To get the slope of the tangent in the (x,y) plane, note that:
    %   x = asph(y)  =>  dx/dy = dAsph(y)
    % so that dy/dx = 1/dAsph(y) (watching for division by zero)
    if abs(dAsp_val) < eps
        slopeTangent = Inf;
    else
        slopeTangent = 1/dAsp_val;
    end
    % The normal is perpendicular to the tangent.
    slopeNormal = -1/slopeTangent;
    
    % Compute the angle of incidence between the ray (with slope "slope")
    % and the normal using:
    %   aInc = atan( -(slopeNormal - slope)/(1 + slopeNormal*slope) )
    aInc = atan( -(slopeNormal - slope) / (1 + slopeNormal * slope) );
end

%% Helper Function: asph
function val = asph(y, r, k, A4, A6, d2)
% asph computes the aspherical lens surface function at y.
%
% The aspherical surface is defined as:
%   x = -y^2/(r*(1 + sqrt(1 - (1+k)*y^2/r^2))) - A4*y^4 - A6*y^6 + d2
%
% If the argument of the square root is negative, the function returns NaN.
    term = 1 - (1+k)*y^2/r^2;
    if term < 0
        val = NaN;
        return;
    end
    s = sqrt(term);
    val = -y^2/(r*(1+s)) - A4*y^4 - A6*y^6 + d2;
end

%% Helper Function: dAsph
function dval = dAsph(y, r, k, A4, A6)
% dAsph computes the derivative of the aspherical surface function with respect to y.
%
% We differentiate:
%   asph(y) = -y^2/(r*(1 + sqrt(1 - (1+k)*y^2/r^2))) - A4*y^4 - A6*y^6 + d2
%
% Common sub-expressions are computed only once. If the square-root argument is negative,
% the function returns NaN.
    alpha = (1+k)/r^2;
    temp = 1 - alpha*y^2;
    if temp < 0
        dval = NaN;
        return;
    end
    s = sqrt(temp);
    % Compute derivative of T(y) = -y^2/(r*(1+s))
    % Let D = r*(1+s). Then T(y) = -y^2/D.
    % Its derivative is computed analytically:
    dT = - (2*y*(1+s) + (alpha*y^3)/s) / (r*(1+s)^2);
    % Derivative of the polynomial terms:
    dPoly = -4*A4*y^3 - 6*A6*y^5;
    dval = dT + dPoly;
end
