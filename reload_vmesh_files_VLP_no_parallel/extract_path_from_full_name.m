function path = extract_path_from_full_name(full_name)



if ismac || isunix
    kkk    = strfind(full_name, '/'); 
elseif ispc
    kkk    = strfind(full_name, '\');  
end

path = full_name(1:kkk(end)-1);









end