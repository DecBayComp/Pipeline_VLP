function path_target_folder = create_a_local_output_folder_and_go_there(name_out)
%%generate folder and go there
%to be found in : Matlab/projet/custom_functions

path      = pwd;
mkdir(name_out);
if ismac||isunix
    full_path =  [path '/' name_out];
elseif ispc
    full_path =  [path '\' name_out]; 
end
cd(full_path);

path_target_folder = full_path;


end