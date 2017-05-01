function [movie_per_frame,assignment]  = add_linkability_number_movie_per_frame(movie_per_frame, name_pipeline, varargin)




try
    warning('off','all');
    pctRunOnAll  warning('off','all'); 
end

if isfield(movie_per_frame, 'nb_linkable')
    assignment = [];
else

    

modal                  = 'hungarian';
[assignment, ~,~ ]     = Assignment_Multi_Mode_Full_Movie(movie_per_frame,modal,[],name_pipeline);
n                      = length(movie_per_frame);

for i = 1 : n-1
    movie_per_frame(i).nb_linkable = length(assignment(i).index_assigned);
    %fprintf('%i\t %i\n', i, length(assignment(i).index_assigned));
end

 movie_per_frame(n) = [];
end













end
