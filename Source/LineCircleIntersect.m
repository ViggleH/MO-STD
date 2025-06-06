%##########################################################################
%                  INTERSECTIONS OF LINE AND CIRCLE
%##########################################################################

%Funtion takes inputs of line(single point on line and slope) and inputs of
%circle(centre point and radius) to determain if and where they intersect.
%Returns the discriminant value( <0 if there is no intersection, 0 if there is a single
%intersection(tangent) and >0 if there is two intersections(secant)). Returns 
%the point(s) of intersection.Returns the more left intersection as (x1,y1) and
%the more right intersetion as (x2,y2)

function [x1,y1,x2,y2,discrim] = LineCircleIntersect(xLine,yLine,m,xCircle,yCircle,radius)

%     a = m^2+1;
%     b = -2*m^2*xLine-2*m*yCircle+2*m*yLine+2*xCircle;
%     c = m^2*xLine^2-2*m*xLine*yCircle-2*m*xLine*yLine-xCircle^2+yCircle^2-2*yCircle*yLine+yLine^2-radius^2; 
 
    g = yLine - m*(xLine); %g is yintercept from line equation: y = mx + b
    a = m^2 + 1;
    b = 2*(m*g - m*yCircle - xCircle);
    c = (yCircle^2 - radius^2 + xCircle^2 - 2*g*yCircle + g^2);

    x1 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
    x2 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
    
    if(x1>x2)
        left=x2;
        x2=x1;
        x1=left;
    end
        
    
    y1 = m*(x1-xLine)+yLine;
    y2 = m*(x2-xLine)+yLine;
    
    discrim = b^2 - 4*a*c;
    
    return
end