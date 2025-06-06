%##########################################################################
%                             SNELL'S LAW
%##########################################################################

%Function takes the index of refraction of incident and transmitted mediums
%with the angle of incidence(in radian) with respect to the normal of the 
%interface and computes the angle of refraction(in radian) with respect to
%the normal of the interface

function [angleRefract,Reflected] = Snells(nIncident,nTransmitted,angleIncident,polAng)
% polAng=45;
    angleRefract = asin((nIncident/nTransmitted)*sin(angleIncident));
    
         %for s polarization component
         Rs_top = (nIncident*cos(angleIncident)-nTransmitted*cos(angleRefract));
         Rs_bot = ((nIncident*cos(angleIncident)+nTransmitted*cos(angleRefract)));
         Rs  = (sin(polAng)^2)*(Rs_top/Rs_bot)^2;
         %for P polarization component
         Rp_top = ((nIncident*cos(angleRefract)-nTransmitted*cos(angleIncident)));
         Rp_bot = ((nIncident*cos(angleRefract)+nTransmitted*cos(angleIncident)));
         Rp  = (cos(polAng)^2)*(Rp_top/Rp_bot)^2;
         
         %total precent lost/reflected
         Reflected =((Rp+Rs)); %if our polAng is 90 isn't that just s-polarization?
         %disp("R = "+Reflected);
    %disp("transmitted = " + (1-Reflected) +" at "+ angleRefract + " radians");
end