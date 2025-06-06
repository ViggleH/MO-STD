function hv3 = computeHypervolume3D(PF, Nadir, Ideal)
    % computeHypervolume3D  Computes the normalized hypervolume for a
    %   three-objective minimization problem using slicing along f1.
    
    % 1) Filter PF so that each point lies within [Ideal, Nadir]
    inBox = all(PF >= Ideal, 2) & all(PF <= Nadir, 2);
    PF    = PF(inBox, :);
    
    if isempty(PF)
        hv3 = 0;
        return;
    end
    
    % 2) Sort PF by the first objective (f1) in ascending order
    PF = sortrows(PF, 1);
    M  = size(PF, 1);
    
    hvAccum    = 0;
    PF2_accum  = zeros(0, 2);  % will collect [f2, f3] up to current slice
    
    for i = 1:M
        % 2a) Add this point’s (f2, f3) to the accumulated set
        PF2_accum(end+1, :) = PF(i, 2:3);  %#ok<AGROW>
        
        % 2b) Compute the raw 2D hypervolume in the (f2, f3) subspace
        %     Reference = [Nadir(2), Nadir(3)], Ideal = [Ideal(2), Ideal(3)]
        area2_i = computeHV2D_unorm(PF2_accum, Nadir(2:3), Ideal(2:3));
        
        % 2c) Determine Δf1 for this slice
        if i < M
            df1 = PF(i+1, 1) - PF(i, 1);
        else
            df1 = Nadir(1) - PF(i, 1);
        end
        
        % 2d) Accumulate slab volume
        hvAccum = hvAccum + area2_i * df1;
    end
    
    % 3) Normalize by the total 3D box volume = ∏_{j=1}^3 (Nadir(j) - Ideal(j))
    totalVol = (Nadir(1) - Ideal(1)) * ...
               (Nadir(2) - Ideal(2)) * ...
               (Nadir(3) - Ideal(3));
    hv3 = hvAccum / totalVol;
end


%% Helper: computeHV2D_unorm
function area2 = computeHV2D_unorm(PF2, N23, I23)
    % computeHV2D_unorm  Computes the raw (unnormalized) hypervolume for a
    %   two-objective minimization problem in the f2–f3 subspace.
    %
    % Steps:
    %   1) Filter PF2 so each point lies within [I23, N23].
    %   2) Extract only the nondominated front (because dominated points
    %      do not contribute to the 2D hypervolume).
    %   3) Sort the nondominated points by f2 ascending.
    %   4) Use the standard “staircase” area formula for 2D HV.
    
    % 1) Filter PF2 into the box [I23, N23]
    inBox2 = all(PF2 >= I23, 2) & all(PF2 <= N23, 2);
    PF2    = PF2(inBox2, :);
    
    if isempty(PF2)
        area2 = 0;
        return;
    end
    
    % 2) Extract the 2D Pareto front (nondominated set) among PF2
    isDominated = false(size(PF2, 1), 1);
    for a = 1:size(PF2,1)
        for b = 1:size(PF2,1)
            if b ~= a
                % If PF2(b,:) dominates PF2(a,:): f2(b) <= f2(a) && f3(b) <= f3(a)
                % (with at least one strict), then mark a as dominated.
                if PF2(b,1) <= PF2(a,1) && PF2(b,2) <= PF2(a,2) && ...
                   (PF2(b,1) < PF2(a,1) || PF2(b,2) < PF2(a,2))
                    isDominated(a) = true;
                    break;
                end
            end
        end
    end
    PF2_nd = PF2(~isDominated, :);
    
    % 3) Sort the nondominated front by f2 ascending
    PF2_nd = sortrows(PF2_nd, 1);
    
    % 4) Compute the “staircase” area under the front to the reference N23
    %    (cf. 2D HV for minimization).
    area2 = (N23(1) - PF2_nd(1,1)) * (N23(2) - PF2_nd(1,2));
    for k = 2:size(PF2_nd, 1)
        width  = N23(1) - PF2_nd(k, 1);
        height = PF2_nd(k-1, 2) - PF2_nd(k, 2);
        area2  = area2 + width * height;
    end
end
