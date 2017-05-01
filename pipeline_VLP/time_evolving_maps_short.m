function time_evolving_maps_short(name_out,movie_per_frame,state,modal, criterion, name_pipeline, fixed_mesh,i, varargin)
 warning('off','all');

if ~exist('fixed_mesh');        fixed_mesh        = []    ; end
if ~exist('name_pipeline');     name_pipeline     = 'main'; end


state = check_state_in_time_evolution(state);
create_a_local_output_folder_and_go_there(name_out);

[movie_per_frame, assignment]  = add_linkability_number_movie_per_frame(movie_per_frame,name_pipeline );
parameters_loc                 = generate_parameters(movie_per_frame, [], name_pipeline);

%% 
[indice_range, n_window] = get_indice_range(criterion,movie_per_frame, name_pipeline);

if nargin == 8
    if ischar(i)
        i = str2num(i);
    end
    
    if i<= n_window
        inum    = num2str(i);
        t_init  = movie_per_frame(indice_range(i,1)).t;
        t_end   = movie_per_frame(indice_range(i,2)).t; 
        nb_tot  = total_number_of_linked_points(movie_per_frame,indice_range(i,1), indice_range(i,2));
        if nb_tot > 3 * parameters_loc.number_per_zone;
            [Maps, Maps_densities]  = Mapping_without_tracking(movie_per_frame(indice_range(i,1):indice_range(i,2)),state,modal,assignment(indice_range(i,1):indice_range(i,2)), fixed_mesh, name_pipeline);
            [Maps]                  = add_features_time_evolving_Maps(Maps,t_init, t_end );
            [Maps]                  = add_one_feature_to_Maps(Maps, nb_tot, 'nb_linked_particles');
        else
            [Maps]           = [];
            [Maps_densities] = [];
        end
        save_individual_Mat_file(Maps, inum);
%         save_individual_density_file(Maps_densities, inum);
        clear Maps;
    end
    
else

    fprintf('n_window %i\n',n_window );
    for i       = 1 : n_window
        inum    = num2str(i);
        t_init  = movie_per_frame(indice_range(i,1)).t;
        t_end   = movie_per_frame(indice_range(i,2)).t; 
        nb_tot  = total_number_of_linked_points(movie_per_frame,indice_range(i,1), indice_range(i,2));
        if nb_tot > 3 * parameters_loc.number_per_zone;
            fprintf('indice range %i\t %i\n',indice_range(i,1), indice_range(i,2) );
            fprintf('state modal %s\t %s\n',state, modal );
            fprintf('name pipeline %s\n',name_pipeline );
            
            
            [Maps, Maps_densities]  = Mapping_without_tracking(movie_per_frame(indice_range(i,1):indice_range(i,2)),state,modal,assignment(indice_range(i,1):indice_range(i,2)), fixed_mesh, name_pipeline);
            [Maps]                  = add_features_time_evolving_Maps(Maps,t_init, t_end );
            [Maps]                  = add_one_feature_to_Maps(Maps, nb_tot, 'nb_linked_particles');
        else
            [Maps]           = [];
            [Maps_densities] = [];
        end
        save_individual_Mat_file(Maps, inum);
%         save_individual_density_file(Maps_densities, inum);
        clear Maps;
    end

end

cd ..;
% kill_parallel_multi_version;



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

