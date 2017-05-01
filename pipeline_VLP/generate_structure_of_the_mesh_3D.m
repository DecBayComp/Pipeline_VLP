function mesh = generate_structure_of_the_mesh_3D(IDX,C,VV,CC,x_tot, y_tot,z_tot, t_tot, dx_tot, dy_tot, dz_tot)




mesh = struct;

for i = 1 : length(C)

    II = IDX == i;
    mesh(i).x          = x_tot(II);
    mesh(i).y          = y_tot(II);
    mesh(i).z          = z_tot(II);
    mesh(i).t          = t_tot(II);
    mesh(i).dx         = dx_tot(II);
    mesh(i).dy         = dy_tot(II);
    mesh(i).dz         = dz_tot(II);

    mesh(i).center_x   = C(i,1);
    mesh(i).center_y   = C(i,2);
    mesh(i).center_z   = C(i,3);

    xx                 = VV(CC{i},1);
    yy                 = VV(CC{i},2);
    zz                 = VV(CC{i},3);


    xx                 = xx(~isinf(xx));
    yy                 = yy(~isinf(yy));
    zz                 = zz(~isinf(zz));


    if (length(xx)<3)

        xx                 = [xx;C(i,1)];
        yy                 = [yy;C(i,2)];  
        zz                 = [yy;C(i,3)];  
        mesh(i).voronoi_x  = xx;
        mesh(i).voronoi_y  = yy;
        mesh(i).voronoi_z  = zz;


    else

        xx                 = [xx;C(i,1)];
        yy                 = [yy;C(i,2)];     
        zz                 = [zz;C(i,3)];     

        K                  = convhull(xx,yy,zz);

        for kkk = 1 : length(K(:,1))

            mesh(i).voronoi_triangular(kkk).x = xx(K(kkk,:)) ;
            mesh(i).voronoi_triangular(kkk).y = yy(K(kkk,:)) ;
            mesh(i).voronoi_triangular(kkk).z = zz(K(kkk,:)) ;

        end


    end

end










end