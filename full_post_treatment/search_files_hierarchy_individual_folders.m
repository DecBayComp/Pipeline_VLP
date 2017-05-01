function [liste, n_liste] = search_files_hierarchy_individual_folders(target_folder,extension, varargin) 


if nargin >= 1
   cd(target_folder); 
end


files   = search_for_files_hierarchy(extension);
n_files = length(files);
liste   = cell(n_files,1);
for i =1 : n_files
    name = files(i).name;
    if ismac||isunix
        kkk    = strfind(name, '/'); 
    elseif ispc
        kkk    = strfind(name, '\');  
    end
    liste{i,1}  = name(1:kkk(end)) ;
end

liste = unique(liste);
n_liste = length(liste);




end