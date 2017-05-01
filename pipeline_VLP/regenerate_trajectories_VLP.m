function trajs = regenerate_trajectories_VLP(movie_per_frame , indice_num, name_pipeline)

[assignment, ~,~ ] = Assignment_Multi_Mode_Full_Movie(movie_per_frame,'hungarian', [], name_pipeline);
[trajs]            = traj_from_assignment(assignment, 2); 
fichier            = fopen(['trajectories_' indice_num '.txt'], 'w');

for l = 1 : length(trajs(:,1))
    fprintf(fichier, '%i\t %f\t %f\t %f\n', trajs(l,1),trajs(l,2),trajs(l,3),trajs(l,4) );
end
fclose(fichier);
            
          


end