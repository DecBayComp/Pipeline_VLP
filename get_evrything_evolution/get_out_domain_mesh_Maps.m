function indice_out = get_out_domain_mesh_Maps(Maps)


n_up          = 10;
n_Maps        = length(Maps);
boolean_z     = isfield(Maps, 'z');
[R_start]     = get_average_point_position_Maps(Maps);
R_centre      = get_centers(Maps, boolean_z, n_Maps);

d_centre      = (R_centre{1} - R_start(1)).^2 + (R_centre{2} - R_start(2) ).^2;
if boolean_z
    d_centre  = d_centre + ( R_centre{3} - R_start(3) ).^2;
end


[~,I]         = sort(d_centre, 'descend');
n_I           = length(I);
I = I(1:min(n_up, n_I));
indice_out = I;


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