function parameters = get_dimension_parameters(movie_per_frame, parameters)

% nb max neighbor dimention
parameters.nb_max_neighbors_2D              = 20;
parameters.nb_max_neighbors_3D              = 40;

%% number of dimentions
if isfield(movie_per_frame, 'z')
    parameters.d                         = 3;
    parameters.nb_max_neighbors          = parameters.nb_max_neighbors_3D;
else
    parameters.d                         = 2;
    parameters.nb_max_neighbors          = parameters.nb_max_neighbors_2D;
end





end