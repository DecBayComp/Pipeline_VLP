function Maps  = embed_value_meighbors_voronoi_3D(Maps)


%   h       = waitbar(0, 'Establish directional neighbors...'); 
nn_Maps = length(Maps);
for i = 1 : nn_Maps

%          waitbar(i/nn_Maps,h);
    nn = length(Maps(i).neighbors);
    x  = Maps(i).center_x;
    y  = Maps(i).center_y;
    z  = Maps(i).center_z;

    indice1 = 1;
    indice2 = 1;
    indice3 = 1;
    indice4 = 1;
    indice5 = 1;
    indice6 = 1;



    %Maps(i).surface = polyarea(Maps(i).tree_x,Maps(i).tree_y);
%                 Maps(i).volume = (Maps(i).border_right - Maps(i).border_left).* (Maps(i).border_up - Maps(i).border_down) .*...
%                                  (Maps(i).border_north - Maps(i).border_south) ;
%               
    tout_shell_loc = [];
    for kk = 1 : length(Maps(i).voronoi_triangular)
         tout_shell_loc = [ tout_shell_loc; Maps(i).voronoi_triangular(kk).x, ...
             Maps(i).voronoi_triangular(kk).y, Maps(i).voronoi_triangular(kk).z];
    end
    [~ ,V_tot]              = convhulln(tout_shell_loc  );
    Maps(i).volume = V_tot;


    if length(Maps(i).x) <=4

        %Maps(i).surface_convhull = 0.;
        Maps(i).volume_convhull = 0.;
    else

        [~ ,V]              = convhulln([Maps(i).x,Maps(i).y,Maps(i).z]  );
        Maps(i).volume_convhull = V;
    end


    Maps(i).plus_z          = [];
    Maps(i).plus_z_index    = [];
    Maps(i).minus_z         = [];
    Maps(i).minus_z_index   = [];
    Maps(i).plus_x          = [];
    Maps(i).plus_x_index    = [];
    Maps(i).minus_x         = [];
    Maps(i).minus_x_index   = [];
    Maps(i).plus_y          = [];
    Maps(i).plus_y_index    = [];
    Maps(i).minus_y         = [];
    Maps(i).minus_y_index   = [];

    for j =1 : nn

        k                               = Maps(i).neighbors(j);
        Maps(i).center_x_neighbors(j,1) = Maps(k).center_x;
        Maps(i).center_y_neighbors(j,1) = Maps(k).center_y;
        Maps(i).center_z_neighbors(j,1) = Maps(k).center_z;

%                     Maps(i).angle_voisin(j,1) =mod( angle(  (Maps(i).center_x_neighbors(j,1) - x) + 1j*(Maps(i).center_y_neighbors(j,1) - y)    ),2*pi);
%                     theta = Maps(i).angle_voisin(j,1);

        [azimuth,elevation,~] = cart2sph(Maps(i).center_x_neighbors(j,1) - x ,...
            Maps(i).center_y_neighbors(j,1) - y , Maps(i).center_z_neighbors(j,1) - z);

        azimuth = mod(azimuth,2*pi);


       % fprintf('azimuth %f\t condition %i\n ',azimuth,  ((azimuth >= 0) && (azimuth <= pi)) );

%                     elseif ( ((azimuth >= 0)       & (azimuth <= 3*pi/8)) |...
%                              ((azimuth >= 13*pi/8) & (azimuth <= 2*pi)) )
%                     elseif ( ((azimuth >= 5*pi/8)  & (azimuth <= 11*pi/8)) )   
%                     elseif ( (azimuth >= pi/8) & (azimuth <= 7*pi/8) )
%                     elseif ( (azimuth >= 9*pi/8) & (azimuth <= 17*pi/8 ) ) 


        if ( elevation >= pi/8 )
            Maps(i).plus_z(indice1)       = Maps(k).center_z;
            Maps(i).plus_z_index(indice1) = k;
            indice1 = indice1 + 1;
        end
        if (elevation <= -pi/8)
            Maps(i).minus_z(indice2)       = Maps(k).center_z;
            Maps(i).minus_z_index(indice2) = k;
            indice2 = indice2 + 1;
        end
        if ( ((azimuth >= 0)       && (azimuth <= pi/2)) ||...
                 ((azimuth >= 3*pi/2) && (azimuth <= 2*pi)) )
             %fprintf('here x %f\n', azimuth);
            Maps(i).plus_x(indice3)       = Maps(k).center_x;
            Maps(i).plus_x_index(indice3) = k;
            indice3 = indice3 + 1;
        end
        if ( ((azimuth >= pi/2)  && (azimuth <= 3*pi/2)) )    
            Maps(i).minus_x(indice4)       = Maps(k).center_x;
            Maps(i).minus_x_index(indice4) = k;
            indice4 = indice4 + 1;
        end
        if ( (azimuth >= 0) && (azimuth <= pi) )
            Maps(i).plus_y(indice5)       = Maps(k).center_y;
            Maps(i).plus_y_index(indice5) = k;
            indice5 = indice5 + 1;
        end
        if ( (azimuth >= pi) && (azimuth <= 2*pi ) ) 
            Maps(i).minus_y(indice6)       = Maps(k).center_y;
            Maps(i).minus_y_index(indice6) = k;
            indice6 = indice6 + 1;
        end

        clear azimuth elevation;

    end

end


total_volume_convhull = 0.;
for i = 1 : nn_Maps
    total_volume_convhull = total_volume_convhull + Maps(i).volume_convhull ;
end

for i = 1 : nn_Maps
    Maps(i).total_volume_convhull = total_volume_convhull;
end

%close(h);



end



