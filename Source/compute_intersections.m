%##########################################################################
%          COMPLETE REFRACTION FUNCTION IS DETECTED ZERO DUPLICATES
%##########################################################################

function [xInts,yInts,intensityProfile] = compute_intersections(geo)

%################## Physical parameters of system #########################
%Can print system paramters for looping

% fprintf('r = %f\n' , r); % number of rays
% fprintf('gel = %f\n' , gel); % refractive index of the gel
% fprintf('span1 = %f\n' , span1); % the span of the fan beam above the horizontal(in radian) (Positive Entry)
% fprintf('span2 = %f\n' , span2); % the span of the fan beam below the horazontal(in radian) (Negative Entry)
% fprintf('edgeDensity = %f\n' , edgeDensity); % Fraction of total ray density incident on the edges of fan beam(.2 radian).
%                                              % Uniform is edgeDesity = (.02+.02)/(span1+(-span2))
% fprintf('dlaser = %f\n' , dlaser); % distance from laser source to acrylic block
% fprintf('doff = %f\n' , doff); % vertical offset of laser source
% fprintf('bath = %f\n' , bath); % refractive index of the matching bath
% fprintf('edgeSpan = %f\n' , edgeSpan); %angle(rad) of increased density on both ends of fan beam

r = geo.r; % number of rays
detN = geo.numDet; % number of detectors

%dlaser = geo.dlaser; % distance from laser source to acrylic block
%doff = geo.doff; % vertical offset of laser source

bath = geo.bath; % refractive index of the matching bath
gel = geo.gel; % refractive index of the gel
air = geo.air;% refractive index of air
acrylic = geo.acrylic; % refractive index of acrylic (block and container)
block = geo.block;

span1 = geo.span1; % the span of the fan beam above the horizontal(in radian) (Positive Entry)
span2 = geo.span2; % the span of the fan beam below the horazontal(in radian) (Negative Entry)

%dia = geo.dia; % 104mm - diameter of bore hole
%wallT = geo.wallT; % 3.5mm - container wall thickness
% hh = geo.hh;

%h = geo.h; 
k = geo.k; % center of bore at (h,k)
h = 0;

%gapT = geo.gapT; % gap thickness i.e. bath thickness
d1 = geo.d1; % distance from front edge of acrylic block to center of bore
d2 = geo.d2; % distance from center of bore hole to detectors, 344mm = total length of acrylic block
x0 = geo.x0; % laser source position on x-axis
det = geo.det; % detectors position
y0 = geo.y0; % laser source position on y-axis
wall = geo.wall; % position of front edge of acrylic block

switch(geo.laserType)
    case 'fan'
        angles = linspace(span1,span2,r); % angle of inclination of each ray from laser source
    case 'optimalCrossover'
        angles = [linspace(-40*pi/180,span2,r/2) linspace(span1,40*pi/180,r/2)]; % angle of inclination of each ray from laser source
    case 'line' %line as opposed to fan
        opp = linspace((d1+d2)*tan(span1),(d1+d2)*tan(span2),r);
        adj = linspace((d1+d2),(d1+d2),length(opp));
        denom = (opp.^2 + adj.^2).^0.5;
        angles = asin(opp./denom);
end

r1 = geo.r1; %gel radius
r2 = geo.r2; %container radius
r3 = geo.r3; %bore radius


% critical angles assuming bath < acrylic
crit1 = asin(bath/block); % critical angle where total internal reflection occurs at Block->Bath interface
crit2 = asin(gel/acrylic); % critical angle where total internal reflection occurs at Acrylic->Gel interface
%crit3 = asin(air/acrylic); % critical angle where total internal reflection occurs at Acrylic->Bath interface
crit4 = asin(air/block); % critical angle where total internal reflection occurs at Block->air interface

%Calculate position of middle of each detector
detSpace = geo.detBayHeight/geo.numDet; %Detectors are 0.8mm appart. Can interpolate detectors
detY = zeros(geo.numDet, 1);
for y = 1:geo.numDet
    detY(y) = detSpace/2 + (geo.detBayHeight/2 - (y*detSpace));
end

%Define block face.  radFace = 0 means flat face
switch(geo.lensType)
    case 'circle'
        hFace = (2*(-d1) + ((2*(-d1))^2 - 4*((-d1)^2 - geo.radFace^2))^0.5)/2;
        kFace = 0;
    case 'ellipse'
        hFace = (-d1)+(geo.bEll^2/(1 + geo.ecc^2))^0.5;
        kFace = 0;
        aEll = (geo.bEll^2/(1 + geo.ecc^2))^0.5;
end

%Define REAR block face.
switch(geo.rearLensType)
    case('conic')
        bEll2 = geo.rearLensGeo(1);
        ecc2 = geo.rearLensGeo(2);
        kFace2 = 0;
        aEll2 = (bEll2^2/(1 + ecc2^2))^0.5;
        hFace2 = d2-aEll2;
    case('poly')
        hFace2 = d2;
        kFace2 = 0;
    case('highPoly')
        hFace2 = d2;
        kFace2 = 0;
    case('aspheric')
        hFace2 = d2;
        kFace2 = 0;
end



%################# Compute Intersection Matrices ##########################

% define size of intersection matrices
% XintersectionMatrix = nan(r,10);
% YintersectionMatrix = nan(r,10);
xInts = zeros(r,10);
yInts = zeros(r,10);
IntMatrix = zeros(r,10);
r0 = 1; %defines the initial ray intesity to be one
rayAng = zeros(r,9); %one less because slopes are inbetween intersection points

% loop over rays
for i = 1:r
    IntMatrix(i,9) = real(r0);
    rayAng(i,9) = angles(i);
    % define laser source
    switch(geo.laserType)
        case 'optimalCrossover'
            xInts(i,10) = x0;
            yInts(i,10) = y0;
        case 'fan'
            xInts(i,10) = x0;
            yInts(i,10) = y0;
        case 'parallel'
            x0 = geo.x0(i);
            y0 = geo.y0(i);
            xInts(i,10) = x0;
            yInts(i,10) = y0;
            angles = zeros(1,r);
        case('crossover')
            x0 = geo.x0(i);
            y0 = geo.y0(i);
            xInts(i,10) = x0;
            yInts(i,10) = y0;
            angles = ones(1,r) * atan(((detY(geo.bottomDetYcross)-(detSpace/2)) - geo.y0(1))/(det - geo.x0(1)));
        case('line')
            xInts(i,10) = x0;
            yInts(i,10) = y0;
    end


    %ray describes the parametertes of the ray [x y slope]
    ray = [x0 y0 tan(angles(i))];

    switch(geo.lensType)
        case 'circle'
            if geo.radFace == 0
                %%%%%% Flat face
                % define intersection points with front edge of acrylic block
                xInts(i,9) = wall;
                yInts(i,9) = ray(3)*(wall - ray(1)) + ray(2);

                % REFRACTION (Air -> AcrylicBlock)
                % call the Snell's Law function to calculate the angle of refraction
                [aRef,iRef] = Snells(air,acrylic,angles(i),geo.polAng);
                rayAng(i,8) = aRef;
                IntMatrix(i,8)=IntMatrix(i,9)*(1-iRef);

                % define parameters of new ray after refraction
                ray = [wall (ray(3)*(wall - ray(1)) + ray(2)) tan(aRef)];

                % call LineCircleIntersect to find intersection points with Acrylic->Bath interface
                [xintl,yintl,~,~,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
                %pick left or right intersection
                xint = xintl;
                yint = yintl;

            else
                %%%%%% Circle face
                % call LineCircleIntersect to find intersection points with Acrylic->Bath interface
                [xintl,yintl,~,~,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),hFace,kFace,geo.radFace);
                %pick left or right intersection
                xint = xintl;
                yint = yintl;

                if discrim > 0

                    % define intersection points with front edge of acrylic block
                    xInts(i,9) = xint;
                    yInts(i,9) = yint;

                    % calculate angle of incidence wrt normal of Air->Acrylic interface
                    aInc = atan((tan(angles(i)) - (yint-kFace)/(xint-hFace))/(1+(yint-kFace)/(xint-hFace)*tan(angles(i))));

                    % REFRACTION (Air -> AcrylicBlock)
                    % call the Snell's Law function to calculate the angle of refraction
                    [p,iRef] = Snells(air,block,aInc,geo.polAng);
                    aRef = p+atan((yint-kFace)/(xint-hFace));
                    rayAng(i,8) = aRef;
                    IntMatrix(i,8)=IntMatrix(i,9)*(1-iRef);
                    ray = [xint yint tan(aRef)];

                    % call LineCircleIntersect to find intersection points with Acrylic->Bath interface
                    [xintl,yintl,~,~,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
                    %pick left or right intersection
                    xint = xintl;
                    yint = yintl;
                end
            end

        case 'ellipse'
            [xintl,yintl,~,~,discrim] = LineEllipseIntersect(ray(1),ray(2),ray(3),hFace,kFace,geo.ecc,geo.bEll);
            %pick left or right intersection
            xint = xintl;
            yint = yintl;

            if discrim > 0

                % define intersection points with front edge of acrylic block
                xInts(i,9) = xint;
                yInts(i,9) = yint;

                % calculate angle of incidence wrt normal of air->acrylic interface
                aInc = atan((tan(angles(i)) - (((aEll^2)*(yint-kFace))/((geo.bEll^2)*(xint-hFace)))) / (1 + (((aEll^2)*(yint-kFace))/((geo.bEll^2)*(xint-hFace)))*tan(angles(i))));

                % REFRACTION (Air -> AcrylicBlock)
                % call the Snell's Law function to calculate the angle of refraction
                [p,iRef] = Snells(air,block,aInc,geo.polAng);
                aRef = -(aInc - angles(i) - p);
                rayAng(i,8) = aRef;
                IntMatrix(i,8)=IntMatrix(i,9)*(1-iRef);
                ray = [xint yint tan(aRef)];

                % call LineCircleIntersect to find intersection points with Acrylic->Bath interface
                [xintl,yintl,~,~,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
                %pick left or right intersection
                xint = xintl;
                yint = yintl;
            end
    end

    % seperate rays that do not intersect Acrylic->Bath interface
    if discrim > 0

        % define intersection points with Acrylic->Bath interface
        xInts(i,8) = xint;
        yInts(i,8) = yint;

        % calculate angle of incidence wrt normal of Acrylic->Bath interface
        aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

        % determain if rays are reflected or refracted
        if abs(aInc) < crit1

            % REFRACTION (AcrylicBlock -> Bath)
            [p,iRef] = Snells(block,bath,aInc,geo.polAng);
            aRef = -(aInc - atan(ray(3)) - p);
            rayAng(i,7) = aRef;
            IntMatrix(i,7)=IntMatrix(i,8)*(1-iRef);
            ray = [xint yint tan(aRef)];
            % calculate intersection points with Bath->Container interface
            [xintl,yintl,~,~,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r2);
            xint = xintl;
            yint = yintl;

            % seperate rays that do not intersect Bath->Container interface
            if discrim > 0

                % define intersection points with Bath->Container Interface
                xInts(i,7) = xint;
                yInts(i,7) = yint;

                % calculate angle of incidence wrt normal of Bath->Container interface
                aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));


                % REFRACTION (Bath -> Container)
                [p,iRef] = Snells(bath,acrylic,aInc,geo.polAng);
                aRef = -(aInc - atan(ray(3)) - p);
                rayAng(i,6) = aRef;
                IntMatrix(i,6)=IntMatrix(i,7)*(1-iRef);
                ray = [xint yint tan(aRef)];
                %calculate intersection points with Container->Gel interface
                [xintl,yintl,~,~,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r1);
                xint = xintl;
                yint = yintl;

                % seperate rays that do not intersect Container->Gel interface
                if discrim > 0

                    % define intersection points with Container->Gel Interface
                    xInts(i,6) = xint;
                    yInts(i,6) = yint;

                    % calculate angle of incidence wrt normal of Container->Gel interface
                    aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

                    if abs(aInc) < crit2

                        % REFRACTION (Container -> Gel)
                        [p,iRef] = Snells(acrylic,gel,aInc,geo.polAng);
                        aRef = -(aInc - atan(ray(3)) - p);
                        rayAng(i,5) = aRef;
                        IntMatrix(i,5)=IntMatrix(i,6)*(1-iRef);
                        ray = [xint yint tan(aRef)];

                        %calculate intersection points with Gel->ContainerBack interface
                        [~,~,xintr,yintr,~] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r1);
                        xint = xintr;
                        yint = yintr;
                        % define intersection points with Gel->ContainerBack
                        xInts(i,5) = xint;
                        yInts(i,5) = yint;

                        % calculate angle of incidence wrt normal of Gel->ContainerBack interface
                        aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

                        % REFRACTION (Gel -> ContainerBack)
                        [p,iRef] = Snells(gel,acrylic,aInc,geo.polAng);
                        aRef = -(aInc - atan(ray(3)) - p);
                        rayAng(i,4) = aRef;
                        IntMatrix(i,4)=IntMatrix(i,5)*(1-iRef);
                        ray = [xint yint tan(aRef)];
                        %calculate intersection points with Container->BathBack interface
                        [~,~,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r2);
                        xint = xintr;
                        yint = yintr;
                        % define intersection points with Container->BathBack
                        xInts(i,4) = xint;
                        yInts(i,4) = yint;

                        % calculate angle of incidence wrt normal of Container->BathBack interface
                        aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

                        if abs(aInc) < crit1

                            % REFRACTION (Container -> BathBack)
                            [p,iRef] = Snells(acrylic,bath,aInc,geo.polAng);
                            aRef = -(aInc - atan(ray(3)) - p);
                            rayAng(i,3) = aRef;
                            IntMatrix(i,3)=IntMatrix(i,4)*(1-iRef);
                            ray = [xint yint tan(aRef)];
                            %calculate intersection points with Bath->AcrylicBack interface
                            [~,~,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
                            xint = xintr;
                            yint = yintr;
                            % define intersection points with Bath->AcrylicBack
                            xInts(i,3) = xint;
                            yInts(i,3) = yint;

                            % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
                            aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

                            % REFRACTION (Bath -> AcrylicBack)
                            [p,iRef] = Snells(bath,block,aInc,geo.polAng);
                            aRef = -(aInc - atan(ray(3)) - p);
                            rayAng(i,2) = aRef;
                            IntMatrix(i,2)=IntMatrix(i,3)*(1-iRef);
                            ray = [xint yint tan(aRef)];
                            %                % define intersection points with Detectors
                            %                XintersectionMatrix(i,2) = det; %supposed to be rear wall
                            %                YintersectionMatrix(i,2) = ray(3)*(det - ray(1)) + ray(2);

                            switch(geo.rearLensType)
                                case('conic')
                                    [~,~,xintr,yintr,discrim] = LineEllipseIntersect(ray(1),ray(2),ray(3),hFace2,kFace2,ecc2,bEll2);
                                    xint = xintr;
                                    yint = yintr;
                                    if discrim > 0
                                        aInc = atan((tan(aRef) - (((aEll2^2)*(yint-kFace2))/((bEll2^2)*(xint-hFace2)))) / (1 + (((aEll2^2)*(yint-kFace2))/((bEll2^2)*(xint-hFace2)))*tan(aRef)));
                                    end
                                case('poly')
                                    [xint,yint,aInc,discrim] = LinePolyIntersect(ray(1),ray(2),ray(3),geo.rearLensGeo);
                                case('highPoly')
                                    [xint,yint,aInc,discrim] = LinePolyIntersect(ray(1),ray(2),ray(3),geo.rearLensGeo);
                                case('aspheric')
                                    [xint,yint, aInc, discrim] = LineAspIntersect(ray(1),ray(2),ray(3),geo.rearLensGeo);
                            end


                            if discrim > 0

                                % define intersection points with back edge of acrylic block
                                xInts(i,2) = xint;
                                yInts(i,2) = yint;

                                if abs(aInc) < crit4 %%%%%%%%%new

                                    % REFRACTION (block -> air)
                                    [p,iRef] = Snells(block,air,aInc,geo.polAng);
                                    aRef = -(aInc - atan(ray(3)) - p);
                                    rayAng(i,1) = aRef;
                                    IntMatrix(i,1)=IntMatrix(i,2)*(1-iRef);
                                    ray = [xint yint tan(aRef)];
                                    %calculate intersection points with Bath->AcrylicBack interface
                                    % define intersection points with Detectors
                                    xInts(i,1) = det; %supposed to be rear wall
                                    yInts(i,1) = ray(3)*(det - ray(1)) + ray(2);

                                else
                                    %Totaly internally reflected rays (Acrylic -> Air)

                                    % REFLECTION (block -> air)
                                    aRef = pi - aInc - (aInc - atan(ray(3)));
                                    ray = [xint yint tan(aRef)];
                                    if aRef > pi
                                        aRef = aRef-2*pi;
                                    end

                                    if aRef > -pi/2 && aRef < pi/2
                                        switch(geo.rearLensType)
                                            case('conic')
                                                [~,~,xintr,yintr,discrim] = LineEllipseIntersect(ray(1),ray(2),ray(3),hFace2,kFace2,geo.ecc2,geo.bEll2);
                                                xint = xintr;
                                                yint = yintr;
                                            case('poly')
                                                [xint,yint,~,discrim] = LinePolyIntersect(ray(1),ray(2),ray(3),geo.rearLensGeo);
                                            case('highPoly')
                                                [xint,yint,~,discrim] = LinePolyIntersect(ray(1),ray(2),ray(3),geo.rearLensGeo);
                                            case('aspheric')
                                                [xint,yint, ~, discrim] = LineAspIntersect(ray(1),ray(2),ray(3),geo.rearLensGeo);
                                        end

                                        if discrim > 0
                                            % define intersection points with back edge of acrylic block
                                            xInts(i,1) = xint;
                                            yInts(i,1) = yint;
                                            rayAng(i,1) = aRef;
                                        end


                                    end
                                end
                            end

                        else
                            %Totaly internally reflected rays (Container -> BathBack)


                        end

                    else
                        %Totaly internally reflected rays (Container -> Gel)

                        % REFLECTION (Container -> Gel)
                        aRef = pi - aInc - (aInc - atan(ray(3)));
                        ray = [xint yint tan(aRef)];

                        %calculate intersection points with Container->BathBack interface
                        [~,~,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r2);
                        xint = xintr;
                        yint = yintr;
                        % define intersection points with Container->BathBack
                        xInts(i,5) = xint;
                        yInts(i,5) = yint;

                        % calculate angle of incidence wrt normal of Container->BathBack interface
                        aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

                        if abs(aInc) < crit1

                            % REFRACTION (Container -> BathBack)
                            [p,iRef] = Snells(acrylic,bath,aInc,geo.polAng);
                            aRef = -(aInc - atan(ray(3)) - p);
                            rayAng(i,4) = aRef;
                            IntMatrix(i,4)=IntMatrix(i,5)*(1-iRef);
                            ray = [xint yint tan(aRef)];
                            %calculate intersection points with Bath->AcrylicBack interface
                            [~,~,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
                            xint = xintr;
                            yint = yintr;
                            % define intersection points with Bath->AcrylicBack
                            xInts(i,4) = xint;
                            yInts(i,4) = yint;

                            % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
                            aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

                            % REFRACTION (Bath -> AcrylicBack)
                            [p,iRef] = Snells(bath,block,aInc,geo.polAng);
                            aRef = -(aInc - atan(ray(3)) - p);
                            rayAng(i,3) = aRef;
                            IntMatrix(i,3)=IntMatrix(i,4)*(1-iRef);
                            ray = [xint yint tan(aRef)];
                            % define intersection points with Detectors
                            xInts(i,3) = det;
                            yInts(i,3) = ray(3)*(det - ray(1)) + ray(2);

                        else
                            %Totaly internally reflected rays (Container -> BathBack)
                            % *Ignore* - will not occur with any geometry
                        end

                    end

                    % rays that go straight through Container
                else

                    % calculate intersection points with Container->BathBack
                    [~,~,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r2);
                    xint = xintr;
                    yint = yintr;
                    % define intersection points of rays that go straight through Container
                    xInts(i,6) = xint;
                    yInts(i,6) = yint;

                    % calculate angle of incidence wrt normal of Container->BathBack interface
                    aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

                    if abs(aInc) < crit1

                        %REFRACTION (Container->BathBack)
                        [p,iRef] = Snells(acrylic,bath,aInc,geo.polAng);
                        aRef = -(aInc - atan(ray(3)) - p);
                        rayAng(i,5) = aRef;
                        IntMatrix(i,5)=IntMatrix(i,6)*(1-iRef);
                        ray = [xint yint tan(aRef)];
                        %calculate intersection points with Bath->AcrylicBack interface
                        [~,~,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
                        xint = xintr;
                        yint = yintr;
                        % define intersection points with Bath->AcrylicBack
                        xInts(i,5) = xint;
                        yInts(i,5) = yint;

                        % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
                        aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

                        %REFRACTION (Bath->AcrylicBack)
                        [p,iRef] = Snells(bath,block,aInc,geo.polAng);
                        aRef = -(aInc - atan(ray(3)) - p);
                        rayAng(i,4) = aRef;
                        IntMatrix(i,4)=IntMatrix(i,5)*(1-iRef);
                        ray = [xint yint tan(aRef)];
                        % define intersection points with Detectors
                        xInts(i,4) = det;
                        yInts(i,4) = ray(3)*(det - ray(1)) + ray(2);


                    else
                        %Totaly internanly reflected rays (Container->BathBack)
                        % Rays that go straight through container wall
                        % *Ignore* - will not occur with any geometry
                    end

                end

                % rays that go straight through Bath
            else

                % calculate intersection points with Bath->AcrylicBack
                [~,~,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
                xint = xintr;
                yint = yintr;
                % define intersection points of rays that go straight through Bath
                xInts(i,7) = xint;
                yInts(i,7) = yint;

                % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
                aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

                % REFRACTION (Bath -> AcrylicBack)
                [p,iRef] = Snells(bath,block,aInc,geo.polAng);
                aRef = -(aInc - atan(ray(3)) - p);
                rayAng(i,6) = aRef;
                IntMatrix(i,6)=IntMatrix(i,7)*(1-iRef);
                ray = [xint yint tan(aRef)];
                % define intersection points with Detectors
                xInts(i,6) = det;
                yInts(i,6) = ray(3)*(det - ray(1)) + ray(2);

            end

        else
            %Totaly internaly reflected rays (AcrylicBlock -> Bath)

            % REFLECTION (AcrylicBlock -> Bath)
            aRef = pi - aInc - (aInc - atan(ray(3)));
            ray = [xint yint tan(aRef)];

            %calculate intersection points with Bath->AcrylicBack interface
            [~,~,xintr,yintr,discrim] = LineCircleIntersect(ray(1),ray(2),ray(3),h,k,r3);
            xint = xintr;
            yint = yintr;
            % define intersection points with Bath->AcrylicBack
            xInts(i,7) = xint;
            yInts(i,7) = yint;

            % calculate angle of incidence wrt normal of Bath->AcrylicBack interface
            aInc = atan((tan(aRef) - (yint-k)/(xint-h))/(1+(yint-k)/(xint-h)*tan(aRef)));

            % REFRACTION (Bath -> AcrylicBack)
            [p,iRef] = Snells(bath,block,aInc,geo.polAng);
            aRef = -(aInc - atan(ray(3)) - p);
            rayAng(i,7) = aRef;
            IntMatrix(i,7)=IntMatrix(i,8)*(1-iRef);
            ray = [xint yint tan(aRef)];
            % define intersection points with Detectors
            xInts(i,6) = det;
            yInts(i,6) = ray(3)*(det - ray(1)) + ray(2);

        end

        % rays that go straight to Detectors
    else

        % define intersection points of rays thats go straight to detectors
        xInts(i,9) = det;
        yInts(i,9) = ray(3)*(det - ray(1)) + ray(2);
        rayAng(i,9) = angles(i);

    end

end

detectorWidth = geo.arrayW; % Width of the detector array
detectorRef = geo.det; % Detector reference
numRays = size(xInts, 1); % Total number of rays

% Initialize intensity profiles
%overallProfile = zeros(1, detN);
intensityProfile = zeros(1, detN);

% Loop through rays and calculate profiles
for rayIndex = 1:numRays
    if nnz(xInts(rayIndex, :)) < 9 % Filter out tangent rays
        continue;
    elseif xInts(rayIndex, 1) == detectorRef
        % Loop through detector bins
        for binIndex = 1:detN
            binUpper = detectorWidth - (binIndex - 1) * (2 * detectorWidth / detN);
            binLower = detectorWidth - binIndex * (2 * detectorWidth / detN);

            % Check if the ray intersects the current bin
            if yInts(rayIndex, 1) <= binUpper && yInts(rayIndex, 1) >= binLower
                %overallProfile(binIndex) = overallProfile(binIndex) + 1;
                intensityProfile(binIndex) = intensityProfile(binIndex) + IntMatrix(rayIndex);
            end
        end
    end
end

end


