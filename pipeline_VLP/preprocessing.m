function [movie_per_frame] = preprocessing(movie_per_frame)
%% ensure that the data is in microns
movie_per_frame = check_units_data(movie_per_frame);



end