function parameters = get_all_non_input_parameters(movie_per_frame, parameters)
format long;

parameters = get_dimension_parameters(movie_per_frame, parameters);
parameters = get_dt_between_frames_parameters(movie_per_frame,parameters);
parameters = get_start_end_time_parameters(movie_per_frame,parameters);
parameters = get_number_movie_parameters(movie_per_frame, parameters);
parameters = get_opt_parameters(parameters);      
parameters = get_param_BP_parameters(parameters)  ; 




end