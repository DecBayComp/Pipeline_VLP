function indice_in = get_mesh_index_inside_sphere_region(Maps, R_start, radius)

n_Maps        = length(Maps);
boolean_z     = isfield(Maps, 'z');
R_centre      = get_centers(Maps, boolean_z, n_Maps);
d_centre      = (R_centre{1} - R_start(1)).^2 + (R_centre{2} - R_start(2) ).^2;
if boolean_z
    d_centre  = d_centre + ( R_centre{3} - R_start(3) ).^2;
end

indice_in = [];
for i = 1 : length(Maps)
   
    if (d_centre(i) < radius^2 )
        indice_in = [indice_in; i];
    end
    
end





end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R_centre = get_centers(Maps, boolean_z,n_Maps)

x_centre          = [];
y_centre          = [];
for i = 1 : n_Maps
    x_centre      = [x_centre; Maps(i).center_x];
    y_centre      = [y_centre; Maps(i).center_y];
end

R_centre{1,1} = x_centre;
R_centre{2,1} = y_centre;


if boolean_z
    z_centre      = [];
    for i = 1 : n_Maps
        z_centre      = [z_centre; Maps(i).center_z];
    end
    R_centre{3,1} = z_centre;

    
end




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



