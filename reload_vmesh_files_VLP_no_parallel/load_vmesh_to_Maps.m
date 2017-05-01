function load_vmesh_to_Maps(full_name_Maps, full_name_vmesh)
warning off;

% full_name_Maps  = [name_Maps '.mat'];
% full_name_vmesh = [name_vmesh '.vmesh'];
% load the files

load(full_name_Maps);
values = dlmread( full_name_vmesh,'\t',9,0 );
    
n_Maps = length(Maps);

x = values(:,3) ;
y = values(:,4) ;
D = values(:,5) ;
V = values(:,8) ;
%V = V - min(V)  ;

for l = 1 : length(x)

    Maps(l).D = D(l,1);
    Maps(l).V = V(l,1);

end
    
[ Maps2 ]           = smooth_maps_multi_methods( Maps, 'diffusion' ,'neighbors');
[ Maps2 ]           = smooth_maps_multi_methods( Maps2, 'potential' ,'neighbors');
[ Maps2 ]           = set_potential_Map_min_to_zeros(Maps2);




Maps = Maps2; 
clearvars -except Maps full_name_Maps;

save(full_name_Maps);

end
