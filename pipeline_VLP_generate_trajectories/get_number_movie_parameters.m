function parameters = get_number_movie_parameters(movie_per_frame, parameters)


%% n tot
parameters.n_movie_per_frame      = length(movie_per_frame) ;
n_tot = 0;
if isfield(movie_per_frame, 'nb_linkable')
    for i = 1 : parameters.n_movie_per_frame
        n_tot = n_tot + movie_per_frame(i).nb_linkable;
    end
    parameters.n_tot = n_tot;
else
    for i = 1 : parameters.n_movie_per_frame
        n_tot = n_tot + length(movie_per_frame(i).x);
    end
    parameters.n_tot = n_tot;
end


end