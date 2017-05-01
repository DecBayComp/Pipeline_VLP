function [Maps] = embed_value_meighbors_tree_3D(Maps)



    
%   h       = waitbar(0, 'Establish directional neighbors...');
nn_Maps = length(Maps);
for i = 1 : nn_Maps

%      waitbar(i/nn_Maps,h);

    nn = length(Maps(i).neighbors);
    x  = Maps(i).center_x;
    y  = Maps(i).center_y;
    z  = Maps(i).center_z;

    indice1 = 1;
    indice2 = 1;
    indice3 = 1;
    indice4 = 1;
    indice5 = 1;
    indice6 = 1;



    Maps(i).volume = (Maps(i).border_right - Maps(i).border_left).* (Maps(i).border_up - Maps(i).border_down) .*...
                     (Maps(i).border_north - Maps(i).border_south) ;


    if length(Maps(i).x) <=4

        Maps(i).volume_convhull = 0.;
    else

       % [~ ,V]              = convhulln(Maps(i).x,Maps(i).y,Maps(i).z  );
        [~ ,V]              = convhulln([Maps(i).x,Maps(i).y,Maps(i).z]  );

        Maps(i).volume_convhull = V;
    end

%             
end


total_volume_convhull = 0.;
for i = 1 : nn_Maps
    total_volume_convhull = total_volume_convhull + Maps(i).volume_convhull ;
end

for i = 1 : nn_Maps
    Maps(i).total_volume_convhull = total_volume_convhull;
end

%close(h);






