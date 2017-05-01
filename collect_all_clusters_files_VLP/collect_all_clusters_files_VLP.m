function collect_all_clusters_files_VLP(target_folder, varargin)

if nargin <1
else
    cd(target_folder);
end

name_folder_cluster      = 'folder_collect_cluster';
files                    = subdir(fullfile(pwd, '*.cluster'));
mkdir(name_folder_cluster);
path                     = pwd;
path_destination         = [path '/' name_folder_cluster];
path_to_collect_clusters = path_destination;

% kkk    = strfind(name, '/');        
% cd(files(i).name(1:kkk(end)-1) );


name_correspondence = struct;
for i = 1 : length(files)
    
    inum = num2str(i);
    name_local = [inum '.cluster'];
    copyfile([files(i).name],[path_destination '/' name_local]);
    
    kk = strfind(files(i).name, '/'); 
    ll = strfind(files(i).name, '.cluster'); 
    
    name_correspondence(i).original_name_and_path = files(i).name;
    name_correspondence(i).original_path = files(i).name(1:kk(end)-1);
    name_correspondence(i).original_name = files(i).name(kk(end)+1:ll(1)-1);
    
    name_correspondence(i).final_name_and_path = [path_destination '/' name_local];
    name_correspondence(i).final_path = path_destination;   
    name_correspondence(i).final_name = inum ;
    
end


cd(path_destination);
save('name_correspondence.mat', 'name_correspondence');
cd(path);

end
