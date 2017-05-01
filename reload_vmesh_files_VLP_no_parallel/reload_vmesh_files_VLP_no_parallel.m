function reload_vmesh_files_VLP_no_parallel(target_folder, varargin)
%% call clear;reload_vmesh_files(1,343,[]);
warning off;


if nargin <1
else
    cd(target_folder);
end

% kill_parallel_multi_version;
% start_check_parallel_multi_version;

files           = subdir(fullfile(pwd, 'Maps*'));
function_load   = @load_vmesh_to_Maps;
function_path   = @extract_path_from_full_name;
function_number = @get_number_from_Maps_files;

% parfor i = 1 : length(files);
for i = 1 : length(files);
    fprintf('%i\n', i);
    try
    name            = files(i).name;
    path            = feval(function_path, name);
    number          = feval(function_number, name);
    number_num      = num2str(number);
    
    if ismac || isunix
        full_name_vmesh = [path '/' number_num '.vmesh'] ;
    elseif ispc
        full_name_vmesh = [path '\' number_num '.vmesh'];
    end

    feval(function_load, name, full_name_vmesh);
    end
end
    %load_vmesh_to_Maps(full_name_Maps, full_name_vmesh);
    %     number          = get_number_from_Maps_files(name);
    %     path            = extract_path_from_full_name(name);