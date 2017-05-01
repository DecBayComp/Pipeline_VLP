function generate_file_for_pure_optimization(Maps, name)
%% generate files from various mesh giving area center positions etc

fichier = fopen([name '.cluster'],'w' );


for i = 1 : length(Maps)
    
    %% area number and positions
    fprintf(fichier, 'ZONE: %i\n',i );
    fprintf(fichier, 'NUMBER_OF_TRANSLOCATIONS: %i\n',length(Maps(i).dx) );
    fprintf(fichier, 'X-CENTRE: %f\n',Maps(i).center_x);
    fprintf(fichier, 'Y-CENTRE: %f\n',Maps(i).center_y);
    fprintf(fichier, 'AREA: %f\n',Maps(i).surface);
    fprintf(fichier, 'AREA_CONVHULL: %f\n',Maps(i).surface_convhull);
    
    
    %% left
    fprintf(fichier, 'NUMBER_OF_LEFT_NEIGHBOURS: %i\n',length(Maps(i).minus_x_index));
    fprintf(fichier, 'LEFT_NEIGHBOURS: ');
    for j = 1 : length(Maps(i).minus_x_index)
       fprintf(fichier, '%i\t',Maps(i).minus_x_index(j)); 
    end
    fprintf(fichier, '\n');
    
    %% right
    fprintf(fichier, 'NUMBER_OF_RIGHT_NEIGHBOURS: %i\n',length(Maps(i).plus_x_index));
    fprintf(fichier, 'RIGHT_NEIGHBOURS: ');
    for j = 1 : length(Maps(i).plus_x_index)
       fprintf(fichier, '%i\t',Maps(i).plus_x_index(j)); 
    end
    fprintf(fichier, '\n');
    
    %% top
    fprintf(fichier, 'NUMBER_OF_TOP_NEIGHBOURS: %i\n',length(Maps(i).plus_y_index));
    fprintf(fichier, 'TOP_NEIGHBOURS: ');
    for j = 1 : length(Maps(i).plus_y_index)
       fprintf(fichier, '%i\t',Maps(i).plus_y_index(j)); 
    end
    fprintf(fichier, '\n');
    
    
    %% bottom
    fprintf(fichier, 'NUMBER_OF_BOTTOM_NEIGHBOURS: %i\n',length(Maps(i).minus_y_index));
    fprintf(fichier, 'BOTTOM_NEIGHBOURS: ');
    for j = 1 : length(Maps(i).minus_y_index)
       fprintf(fichier, '%i\t',Maps(i).minus_y_index(j)); 
    end
    fprintf(fichier, '\n');
    
    %% displacement
    fprintf(fichier, 'DX\t DY\t DT\n');
    for j = 1 : length(Maps(i).dx)
        fprintf(fichier, '%f\t %f\t %f\n', Maps(i).dx(j),Maps(i).dy(j),Maps(i).parameters.dt_theo );
    end
    fprintf(fichier, '\n\n');
    
    
    
end


%% The End
fclose(fichier);

    
    





