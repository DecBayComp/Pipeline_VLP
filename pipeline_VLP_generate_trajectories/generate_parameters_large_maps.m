function parameters = generate_parameters_large_maps(movie_per_frame, parameters, varargin)


parameters = get_all_non_input_parameters(movie_per_frame, parameters);

    parameters = get_tree_mesh_parameters(parameters);
    %% positioning noise
    parameters.sigma                    = 0.025;

    %% diffusion
    parameters.D_high                   = 0.3;
    parameters.length_high              = 2.*sqrt(2*parameters.d*parameters.D_high*parameters.dt_theo);

    %% mesh information
    parameters.number_per_zone          = 50;
    parameters.number_of_cluster        = floor( parameters.n_tot./parameters.number_per_zone )+1;

    %% boosting
    parameters.number_per_zone_mean     = 50;%useless for you
    parameters.number_per_zone_sigma    = 10;%useless for you
    parameters.number_per_zone_min      = 15;
    parameters.number_per_zone_max      = 40;
    % resolution of the final map
    parameters.Maps_dx                  = 0.15;
    % nb mesh generated 
    parameters.nb_trial                 = 5;
    
    % density estimation
    parameters.density_dx               = parameters.Maps_dx ;
    parameters.desnity_sigma            = 0.30;
    parameters.density_sigma_noise      = 0.015;
    %parameters                          = get_Gaussian_desnity_parameters(parameters);
    %% time evolution parameters
    parameters.duration                 = 30    ;
    parameters.t_sliding                = 10  ;
    parameters.nb_points                = 100000;
    parameters.nb_sliding               = 50000;
    %% nb of processors
    parameters.processors               = 6;
    
    
end



















