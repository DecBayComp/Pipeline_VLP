function [ mesh ] = voronoi_meshing_3D( tout, parameters,initial_mesh, varargin  )
%function [ mesh ,VV,CC,C] = voronoi_meshing_3D( tout, parameters,initial_mesh, varargin  )


    
[x_tot, y_tot, z_tot,  t_tot, dx_tot , dy_tot, dz_tot]  = build_point_map(tout);
%%
if isempty(initial_mesh)
   nargin = 2; 
end    



if nargin == 2

    number_of_cluster                       = parameters.number_of_cluster;
    R                                       = [x_tot, y_tot,z_tot];    
    [initial_C]                             = initialize_voronoi_bubble(R, number_of_cluster);
    number_of_cluster                       = length(initial_C(:,1));
    parameters.number_of_cluster            = number_of_cluster;
    [IDX, C, ~]                             = kmeans(R, number_of_cluster, 'Start',initial_C, 'EmptyAction', 'singleton','OnlinePhase','off','Options', parameters.opts);
    nb_max_neighbors                        = parameters.nb_max_neighbors_3D;
    [VV,CC]                                 = voronoin(C) ;

    mesh                                    = generate_structure_of_the_mesh_3D(IDX,C,VV,CC,x_tot, y_tot,z_tot, t_tot, dx_tot, dy_tot, dz_tot);
    mesh                                    = give_neighbors(mesh, C , nb_max_neighbors );

elseif nargin ==3

    mesh = initial_mesh;
    for i = 1 : length(initial_mesh)

        II                 = inpolygon(x_tot,y_tot,initial_mesh(i).voronoi_x,initial_mesh(i).voronoi_y);
        mesh(i).x          = x_tot(II);
        mesh(i).y          = y_tot(II);
        mesh(i).z          = z_tot(II);

        mesh(i).t          = t_tot(II);
        mesh(i).dx         = dx_tot(II);
        mesh(i).dy         = dy_tot(II);
        mesh(i).dz         = dz_tot(II);



    end




end
        
        
       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
