function parameters = generate_parameters_neurons(movie_per_frame, parameters, varargin)


parameters = get_all_non_input_parameters(movie_per_frame, parameters);

    parameters = get_tree_mesh_parameters(parameters);
    %% positioning noise
    parameters.sigma                    = 0.025;

    %% diffusion
    parameters.D_high                   = 0.25;
    parameters.length_high              = 2.*sqrt(2*parameters.d*parameters.D_high*parameters.dt_theo);

    %% mesh information
    parameters.number_per_zone          = 25;
    parameters.number_of_cluster        = floor( parameters.n_tot./parameters.number_per_zone )+1;
    if parameters.number_of_cluster >200
        number_of_cluster = 200;
        parameters.number_per_zone     = floor(parameters.n_tot./number_of_cluster);
        parameters.number_of_cluster   = floor( parameters.n_tot./parameters.number_per_zone )+1;
    end
    
    parameters.min_number_per_zone      = 8;
    parameters.boolean_min              = 1;
%     parameters.mode_initialisation_mesh = '';
    
    %% boosting
    parameters.number_per_zone_mean     = 50;%useless for you
    parameters.number_per_zone_sigma    = 10;%useless for you
    parameters.number_per_zone_min      = 15;
    parameters.number_per_zone_max      = 40;
    % resolution of the final map
    parameters.Maps_dx                  = 0.01;
    % nb mesh generated 
    parameters.nb_trial                 = 5;
    
    % density estimation
    parameters.density_dx               = parameters.Maps_dx ;
    parameters.desnity_sigma            = 0.15;
    parameters.density_sigma_noise      = 0.1;
    %parameters                          = get_Gaussian_desnity_parameters(parameters);
    %% time evolution parameters
    parameters.duration                 = 1    ;
    parameters.t_sliding                = 0.5  ;
    parameters.nb_points                = 1000;
    parameters.nb_sliding               = 100;
    %% nb of processors
    parameters.processors               = 6;
    
    
    
    
end



















