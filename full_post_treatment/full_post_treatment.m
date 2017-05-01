function full_post_treatment(i, target_folder,lambda,value, varargin)


if ischar(i)
     i = str2num(i);
end
    
%% parameters for analysis
% value            = 2e4;
crit_density     = value;
factor           = 2;
dt_theo          = 20e-3;
sigma            = 30e-3;


%% listing files 
path             = pwd;
copy_sub_trajs_in_folders(path);
extension        = 'vmesh';
[liste, n_liste] = search_files_hierarchy_individual_folders(target_folder,extension) ;

if ( i<=n_liste )
% for i = 1 : n_liste
    cd(liste{i});
    [indice_min, indice_max] = get_limit_index_Maps();
    files_loc = dir('trajectories*');
    trajs = load(files_loc(1).name);
    [x_start, y_start, densities, XX, YY, radius] = set_seed_densities_multi_VLP(trajs, value);
    n_start = length(x_start);
    for j = 1 : n_start
        fprintf('%i\t %i\n', j, n_start);
        x_centre_cluster = x_start(j);
        y_centre_cluster = y_start(j);
        radius_loc      = radius(j);
        [t, density]  = get_evolution_density_in_VLP(x_start(j), y_start(j), indice_min, indice_max, radius(j), factor);
        [t, diffusivity,  diffusivity_around, diffusivity_variation]       = get_evolution_diffusivity_VLP(x_start(j), y_start(j), indice_min, indice_max, radius(j), factor);
        [t, potential_depth,  potential_depth_far, potential_depth_global] = get_evolution_potential_depth_VLP(x_start(j), y_start(j), indice_min, indice_max, radius(j), factor, lambda, sigma, dt_theo);
        [t, free_energy,  free_energy_all, free_energy_variation]          = get_evolution_free_energy_VLP(x_start(j), y_start(j), indice_min, indice_max, radius(j), factor, lambda, sigma, dt_theo);
        [t, entropy,  entropy_variation, entropy_all]                      = get_evolution_entropy_VLP(x_start(j), y_start(j), indice_min, indice_max, radius(j), factor, lambda, sigma, dt_theo);
        j_num = num2str(j);
        save(['cluster_' j_num '.mat'],'t', 'density', 'diffusivity',  'diffusivity_around', 'diffusivity_variation', ...
         'potential_depth',  'potential_depth_far', 'potential_depth_global', ...
        'free_energy',  'free_energy_all', 'free_energy_variation',   'entropy',  'entropy_variation', 'entropy_all', ...
        'trajs', 'x_centre_cluster', 'y_centre_cluster', 'radius_loc', 'factor','crit_density', 'indice_min', 'indice_max',...
        'XX', 'YY', 'densities');
        clear('t' , 'density', 'diffusivity',  'diffusivity_around', 'diffusivity_variation', ...
         'potential_depth',  'potential_depth_far', 'potential_depth_global', ...
        'free_energy',  'free_energy_all', 'free_energy_variation',   'entropy',  'entropy_variation', 'entropy_all', ...
         'x_centre_cluster', 'y_centre_cluster');
       
    end
    
    
end







cd(path);


end