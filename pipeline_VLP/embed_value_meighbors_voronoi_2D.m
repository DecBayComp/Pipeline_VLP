function Maps = embed_value_meighbors_voronoi_2D(Maps)

%fprintf('here\n');

%   h       = waitbar(0, 'Establish directional neighbors...');
nn_Maps = length(Maps);
if ~isfield(Maps,'neighbors');

else

    for i = 1 : nn_Maps
    %      waitbar(i/nn_Maps,h);


        nn = length(Maps(i).neighbors);
        x  = Maps(i).center_x;
        y  = Maps(i).center_y;

        indice1 = 1;
        indice2 = 1;
        indice3 = 1;
        indice4 = 1;

        % fprintf('%i\t %i\n ',i, length(Maps(i).x));

        Maps(i).surface = polyarea(Maps(i).voronoi_x,Maps(i).voronoi_y);
        if length(Maps(i).x) <=3

            Maps(i).surface_convhull = 0.;
        else
            try
                K               = convhull(Maps(i).x,Maps(i).y );
                Maps(i).surface_convhull = polyarea(Maps(i).x(K),Maps(i).y(K) );
            catch 
                 Maps(i).surface_convhull = 0.; 
            end
        end

        for j =1 : nn

            k = Maps(i).neighbors(j);
            Maps(i).center_x_neighbors(j,1) = Maps(k).center_x;
            Maps(i).center_y_neighbors(j,1) = Maps(k).center_y;
            Maps(i).angle_voisin(j,1) =mod( angle(  (Maps(i).center_x_neighbors(j,1) - x) + 1j*(Maps(i).center_y_neighbors(j,1) - y)    ),2*pi);
            theta = Maps(i).angle_voisin(j,1);


            if  (( theta>7*pi/4 )||( theta<pi/4 ))
                Maps(i).plus_x(indice1)       = Maps(k).center_x;
                Maps(i).plus_x_index(indice1) = k;
                indice1 = indice1 + 1;
            elseif  (( theta>pi/4 )&&( theta<3*pi/4 ))
                Maps(i).plus_y(indice2)       = Maps(k).center_y;
                Maps(i).plus_y_index(indice2) = k;
                indice2 = indice2 + 1;
            elseif  (( theta>3*pi/4 )&&( theta<5*pi/4 ))
                Maps(i).minus_x(indice3)       = Maps(k).center_x;
                Maps(i).minus_x_index(indice3) = k;
                indice3 = indice3 + 1;
            elseif (( theta>5*pi/4 )&&( theta<7*pi/4 ))
                Maps(i).minus_y(indice4)       = Maps(k).center_y;
                Maps(i).minus_y_index(indice4) = k;
                indice4 = indice4 + 1;
            end

        end

    end

    total_surface_convhull = 0.;
    for i = 1 : nn_Maps
        total_surface_convhull = total_surface_convhull + Maps(i).surface_convhull ;
    end

    for i = 1 : nn_Maps
        Maps(i).total_surface_convhull = total_surface_convhull;
    end 
end



%close(h);













end