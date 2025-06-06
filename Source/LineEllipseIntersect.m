%##########################################################################
%                  INTERSECTIONS OF LINE AND ELLIPSE
%##########################################################################

%Funtion takes inputs of line(single point on line and slope) and inputs of
%ellipse(centre point, eccentricity and semi-major axis) to determain if and where they intersect.
%Returns the discriminant value( <0 if there is no intersection, 0 if there is a single
%intersection(tangent) and >0 if there is two intersections(secant)). Returns 
%the point(s) of intersection.Returns the more left intersection as (x1,y1) and
%the more right intersetion as (x2,y2), that is, x1<x2.

function [x1,y1,x2,y2,discrim] = LineEllipseIntersect(xLine,yLine,slope,xEllipse,yEllipse,eccEllipse,bEllipse)

  aEllipse = (bEllipse^2/(1 + eccEllipse^2))^0.5;


  a = 1/aEllipse^2 + slope^2/bEllipse^2;
  b = -2*xEllipse/aEllipse^2+(2*(-slope*xLine+yLine-yEllipse))*slope/bEllipse^2;
  c = (xEllipse^2/aEllipse^2+(-slope*xLine+yLine-yEllipse)^2/bEllipse^2) - 1; 

  x1 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
  x2 = (-b + sqrt(b^2 - 4*a*c))/(2*a);  

  if(x1>x2)
      left=x2;
      x2=x1;
      x1=left;
  end

  y1 = slope*(x1-xLine)+yLine;
  y2 = slope*(x2-xLine)+yLine;

  discrim = b^2 - 4*a*c;


  return
  
end