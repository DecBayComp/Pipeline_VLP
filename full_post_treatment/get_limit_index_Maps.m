function [indice_min, indice_max] = get_limit_index_Maps()


files_loc = dir('Maps*');
index = [];
for i = 1 : length(files_loc)
   number      = get_number_from_trajectories_files(files_loc(i).name); 
    index = [index; number];
    
end


indice_min = min(index);
indice_max = max(index);



end