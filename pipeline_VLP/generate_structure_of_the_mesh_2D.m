function mesh = generate_structure_of_the_mesh_2D(IDX,C,VV,CC,x_tot, y_tot, t_tot, dx_tot, dy_tot)


mesh = struct;

for i = 1 : length(C)

    II                 = IDX == i;
    mesh(i).x          = x_tot(II);
    mesh(i).y          = y_tot(II);
    mesh(i).t          = t_tot(II);
    mesh(i).dx         = dx_tot(II);
    mesh(i).dy         = dy_tot(II);
    mesh(i).center_x   = C(i,1);
    mesh(i).center_y   = C(i,2); 
    xx                 = VV(CC{i},1);
    yy                 = VV(CC{i},2);
    xx                 = xx(~isinf(xx));
    yy                 = yy(~isinf(yy));

    if (length(xx)<3)

        xx                 = [xx;C(i,1)];
        yy                 = [yy;C(i,2)];      
        mesh(i).voronoi_x  = xx;
        mesh(i).voronoi_y  = yy;

    else

        xx                 = [xx;C(i,1)];
        yy                 = [yy;C(i,2)];     
        K                  = convhull(xx,yy);
        mesh(i).voronoi_x  = xx(K);
        mesh(i).voronoi_y  = yy(K);

    end

end










end