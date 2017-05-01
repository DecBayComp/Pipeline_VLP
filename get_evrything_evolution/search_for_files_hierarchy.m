function files = search_for_files_hierarchy(extension)


files          = subdir(fullfile(pwd, ['*.' extension]));



end