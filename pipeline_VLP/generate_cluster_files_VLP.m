function generate_cluster_files_VLP(name_folder)

function_number = @get_number_from_Maps_files;
path            = pwd;
target          = [path '/' name_folder];

% try
cd (target);
files = dir('Maps*');

for j = 1 : length(files)
    jnum = num2str(j);
    name = files(j).name;
    load(name);
    number          = feval(function_number, name);
    number_num      = num2str(number);
    try
    generate_file_for_pure_optimization(Maps, number_num);
    catch
        
    end
    clear Maps;
end
cd (path);
% catch
%     
% end

        






















