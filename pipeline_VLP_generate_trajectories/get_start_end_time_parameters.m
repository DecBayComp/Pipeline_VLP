function parameters = get_start_end_time_parameters(movie_per_frame,parameters)

% start end
parameters.t_init   = movie_per_frame(1).t;
parameters.t_end    = movie_per_frame(end).t;

end