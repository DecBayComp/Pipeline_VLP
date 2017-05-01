function destroy_files_evrywhere_hierarchy(target_folder, name)
 
warning off;

cd(target_folder);

files_to_destroy  =  subdir(fullfile(pwd, name));
for i = 1 : length(files_to_destroy);
    try delete(files_to_destroy(i).name); end
    try system(['rm -r ' files_to_destroy(i).name]); end
end;


end