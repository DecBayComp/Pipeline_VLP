function clusters = cluster_to_movie_per_frame_VLP(clusters, t_min, t_max, dt_theo)

for i = 1 : length(clusters)
    movie_per_frame = convert_traj_struct(clusters(i).xyt); 
    movie_per_frame = regularize_movie_per_frame_clusters_VLP(movie_per_frame, t_min, t_max, dt_theo);
    clusters(i).movie_per_frame = movie_per_frame ;
end







end