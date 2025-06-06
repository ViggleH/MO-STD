function [xint, yint, aInc, discrim] = LineHighPolyIntersect(xIntercept, yIntercept, slope, polyCoeffs)
% LinePolyIntersect finds the intersection of a ray with a lens defined by a polynomial.
%
%   [xint, yint, aInc, discrim] = LinePolyIntersect(xIntercept, yIntercept, slope, polyCoeffs)
%
% Inputs:
%   xIntercept, yIntercept, slope : define the ray as 
%         y = slope*(x - xIntercept) + yIntercept
%   polyCoeffs : coefficients of the polynomial that defines the lens
%                (so that x = polyval(polyCoeffs, y))
%
% Outputs:
%   xint, yint : coordinates of the intersection point (if found)
%   aInc       : angle of incidence (radians) between the ray and the lens normal
%   discrim    : indicator (1 if intersection is valid, -1 otherwise)
%
% Improved Performance: Instead of using fzero, this version forms a polynomial
% equation in y and uses the roots function to find the intersection(s).

    % Determine the degree (n) of the lens polynomial.
    n = length(polyCoeffs) - 1;
    
    % --- Form the polynomial for the ray ---
    % The ray is given by:
    %   x = (y - yIntercept)/slope + xIntercept.
    % Expressed as a polynomial in y, this is a linear function:
    %   (1/slope)*y + (xIntercept - yIntercept/slope).
    % We create a vector of coefficients of length (n+1) by padding with zeros:
    polyRay = [zeros(1, n-1), 1/slope, (xIntercept - yIntercept/slope)];
    
    % --- Form the intersection polynomial ---
    % We require: polyval(polyCoeffs, y) - polyval(polyRay, y) = 0.
    polyDiff = polyCoeffs - polyRay;
    
    % Solve for all roots of the polynomial equation.
    r = roots(polyDiff);
    
    % Keep only the real roots.
    realRoots = r(imag(r)==0);
    
    % If no real intersection exists, return NaN values.
    if isempty(realRoots)
        xint = NaN;
        yint = NaN;
        discrim = -1;
        aInc = NaN;
        return;
    end
    
    % Choose the root that is closest to the ray's yIntercept.
    [~, idx] = min(abs(realRoots - yIntercept));
    yint = realRoots(idx);
    
    % Evaluate the lens polynomial at yint to get the intersection x-coordinate.
    xint = polyval(polyCoeffs, yint);
    
    % Use a simple test on xint (this may be application‐specific).
    if xint > 0
        discrim = 1;
    else
        xint = NaN;
        discrim = -1;
    end
    
    % --- Compute the normal and the incidence angle ---
    % The lens is defined as x = poly(y), so the derivative dx/dy is:
    dpoly = polyder(polyCoeffs);
    dxdY = polyval(dpoly, yint);
    
    % The slope of the tangent in the (x,y) plane is the reciprocal: dy/dx = 1/(dx/dy)
    % (Note: Check for division by zero if needed.)
    slopeTangent = 1 / dxdY;
    
    % The normal is perpendicular to the tangent.
    slopeNormal = -1 / slopeTangent;
    
    % Compute the angle of incidence between the ray (slope) and the normal.
    % The angle between two lines with slopes m1 and m2 is given by:
    %   theta = atan(abs((m2 - m1)/(1 + m1*m2))).
    % Here, we use the sign convention as in your original code.
    aInc = atan( -(slopeNormal - slope) / (1 + slopeNormal*slope) );
    
end
