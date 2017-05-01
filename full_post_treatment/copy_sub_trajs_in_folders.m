function copy_sub_trajs_in_folders(target_folder, varargin)


cd(target_folder);

files             = subdir(fullfile(pwd, 'trajectories_*'));
n_files           = length(files);


for i  = 1 : n_files
    
    try
        if ismac||isunix
            kkk           = strfind(files(i).name, '/'); 
        elseif ispc
            kkk           = strfind(files(i).name, '\');  
        end

        name_folder = files(i).name(1:kkk(end)-1);
        cd(name_folder);
        number      = get_number_from_trajectories_files(files(i).name);
        number_num  = num2str(number);
        name_file   = files(i).name(kkk(end)+1:end);

        original_path = name_folder;
        if ismac||isunix
            final_path    = [name_folder '/' number_num];
        elseif ispc
            final_path    = [name_folder '\' number_num];
        end
        original_name = name_file;
        final_name    = name_file;
        copy_something_somewhere_change_name(original_path, final_path, original_name, final_name)
    end
    
end

cd(target_folder);


end