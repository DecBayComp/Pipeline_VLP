function number = get_number_from_trajectories_files(name)

try
if ismac || isunix
    kkk    = strfind(name, '/'); 
elseif ispc
    kkk    = strfind(name, '\');  
end

name   = name(kkk(end)+1:end);
end
kkk    = strfind(name, '.');
lll    = strfind(name, '_');
number = str2double(name(lll+1:kkk-1));






end







