function spread_vmesh_to_origin_files_neurons_no_parallel(target_folder, varargin)

warning off;

cd(target_folder);


load('name_correspondence.mat');
n_name_correspondance = length(name_correspondence);
% files = dir('*.cluster');
path  = pwd;
function_copy = @copy_something_somewhere_change_name;
% 
% kill_parallel_multi_version;
% start_check_parallel_multi_version;

original_name = cell(n_name_correspondance,1);
final_name    = cell(n_name_correspondance,1);
original_path = cell(n_name_correspondance,1);
for i = 1: n_name_correspondance
    %fprintf('%i\n', i);
    original_name{i} = name_correspondence(i).original_name;
    final_name{i}      = name_correspondence(i).final_name;
    original_path{i}   = name_correspondence(i).original_path;
end


for i = 1 : length(name_correspondence)
    try
    fprintf('%i\n', i);
    original_rename = [original_name{i} '.vmesh']; 
    final_rename    = [final_name{i} '.vmesh'];
    feval(function_copy,path,original_path{i}, final_rename, original_rename);
    end
end





end