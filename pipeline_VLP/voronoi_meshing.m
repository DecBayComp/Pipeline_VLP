function [ mesh ] = voronoi_meshing( tout, parameters,initial_mesh, varargin  )
%% Voronoi Meshing


number_var = nargin;
[x_tot, y_tot, t_tot,dx_tot , dy_tot]   = build_point_map(tout);
%%
if isempty(initial_mesh)
   number_var = 2; 
end

if number_var == 2

    number_of_cluster                       = parameters.number_of_cluster;
    boolean_min                             = parameters.boolean_min ;
    min_number_per_zone                     = parameters.min_number_per_zone;
    R                                       = [x_tot, y_tot];
    [initial_C]                             = initialize_voronoi_bubble(R, number_of_cluster);    
%     [initial_C]                             = initialize_voronoi_random(R, number_of_cluster);
    number_of_cluster                       = length(initial_C(:,1));
    parameters.number_of_cluster            = number_of_cluster;
%     fprintf('number of cluster %i\t number of points %i\n', number_of_cluster, length(R(:,1)) );
    if boolean_min
        [IDX, C,parameters]                 = corrected_kmeans_fuse_outliers_area(R, initial_C, parameters);
    else
        [IDX, C, ~]                         = kmeans(R, number_of_cluster, 'Start',initial_C, 'EmptyAction', 'singleton','OnlinePhase','off','Options', parameters.opts);
    end
    
    
    nb_max_neighbors                        = parameters.nb_max_neighbors_2D;
    [VV,CC]                                 = voronoin(C) ;

    mesh = generate_structure_of_the_mesh_2D(IDX,C,VV,CC,x_tot, y_tot, t_tot, dx_tot, dy_tot);
    [mesh ]                                 = give_neighbors(mesh, C,nb_max_neighbors );

elseif number_var ==3

    mesh = initial_mesh;
    for i = 1 : length(initial_mesh)

        II                 = inpolygon(x_tot,y_tot,initial_mesh(i).voronoi_x,initial_mesh(i).voronoi_y);
        mesh(i).x          = x_tot(II);
        mesh(i).y          = y_tot(II);
        mesh(i).t          = t_tot(II);
        mesh(i).dx         = dx_tot(II);
        mesh(i).dy         = dy_tot(II);


    end




end
        
        
       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
