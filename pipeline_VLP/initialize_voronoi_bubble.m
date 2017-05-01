function initial_C =  initialize_voronoi_bubble(R, number_of_cluster)
%% initialize the position of the center points of voronoi meshing
% with the bubbling scheme
%to be found in : Matlab/projet/Mapping_without_tracking/voronoi



[~, m_R] = size(R);

if m_R ==2
    
    x_tot         = R(:,1);
    y_tot         = R(:,2);
    inc           = 1e-3;
    [x_min ,I]    = min(x_tot,[],1);
    y_min         = y_tot(I);
    totalCount    = 0;
    n_tot         = length(x_tot);
%     pointsPerZone = round(n_tot /number_of_cluster );
    pointsPerZone = ceil(n_tot /number_of_cluster );
    initial_C     = zeros(number_of_cluster,2);
    indice        = 1;

    while (n_tot>0)

        zoneCount   = 0;
        radius      = 0.0;

        while( (zoneCount < pointsPerZone)  )

            radius      = radius + inc;
            II          = (x_tot - x_min).^2 + (y_tot - y_min).^2 < radius.^2; 
            zoneCount   = sum(II);

            if ( (zoneCount >= pointsPerZone) || (zoneCount == n_tot) )

                if (zoneCount >pointsPerZone ) 
                   JJ = find(II == 1);
                   II(JJ(pointsPerZone+1:zoneCount)) = 0; 
                end

                initial_C(indice,1) = mean(x_tot(II));
                initial_C(indice,2) = mean(y_tot(II));       
                x_tot               = x_tot(~II);
                y_tot               = y_tot(~II);
                [x_min ,I]          = min(x_tot,[],1);
                y_min               = y_tot(I);
                n_tot               = length(x_tot);
                indice              = indice + 1;

                break;

            end


        end

        totalCount = totalCount + zoneCount;

    end

    
    
elseif m_R == 3
      
    
    
    x_tot         = R(:,1);
    y_tot         = R(:,2);
    z_tot         = R(:,3);
    inc           = 1e-3;
    [x_min ,I]    = min(x_tot,[],1);
    y_min         = y_tot(I);
    z_min         = z_tot(I);

    totalCount    = 0;
    n_tot         = length(x_tot);
    pointsPerZone = round(n_tot /number_of_cluster );
    initial_C     = zeros(number_of_cluster,2);
    indice = 1;

    % h = waitbar(0, 'initializing voronoi clusters...');

    while (n_tot>0)

        zoneCount   = 0;
        radius      = 0.0;

        while( (zoneCount < pointsPerZone)  )

            radius      = radius + inc;
            II          = (x_tot - x_min).^2 + (y_tot - y_min).^2 + (z_tot - z_min).^2< radius.^2; 
            zoneCount   = sum(II);

            if ( (zoneCount >= pointsPerZone) || (zoneCount == n_tot) )

                if (zoneCount >pointsPerZone ) 
                   JJ = find(II == 1);
                   II(JJ(pointsPerZone+1:zoneCount)) = 0; 
                end


                initial_C(indice,1) = mean(x_tot(II));
                initial_C(indice,2) = mean(y_tot(II));       
                initial_C(indice,3) = mean(z_tot(II));       

                x_tot               = x_tot(~II);
                y_tot               = y_tot(~II);
                z_tot               = z_tot(~II);

                [x_min ,I]          = min(x_tot,[],1);
                y_min               = y_tot(I);
                z_min               = z_tot(I);

                n_tot               = length(x_tot);
                indice              = indice + 1;

                break;

            end


        end

        totalCount = totalCount + zoneCount;

    end



    
    
end


end