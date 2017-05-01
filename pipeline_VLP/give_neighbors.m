function [mesh ] = give_neighbors(mesh, C,nb_max_neighbors )
%% give the meighbors in a voronoi mesh
%to be found in : Matlab/projet/Mapping_without_tracking/voronoi

[~,m_C] = size(C);

if m_C ==2 

    x               = C(:,1);
    y               = C(:,2);
    [d]             = give_distance_matrix(x, y);
    [~,I]           = sort(d,2,'ascend');
    n_I             = length(I(1,:));
    n_eff_neighbors = min(n_I-1,nb_max_neighbors);
    I               = I(:, 2 : n_eff_neighbors + 1 );
    n_mesh          = length(mesh);
    
    for i =1: n_mesh
        
        indice = 1;
        R_i  = [mesh(i).voronoi_x(1:end-1), mesh(i).voronoi_y(1:end-1),mesh(i).voronoi_x(2:end), mesh(i).voronoi_y(2:end) ] ;
        
        for j = 1 : n_eff_neighbors
           
            R_j_1  = [mesh(I(i,j)).voronoi_x(1:end-1)  , mesh(I(i,j)).voronoi_y(1:end-1)  , mesh(I(i,j)).voronoi_x(2:end)      , mesh(I(i,j)).voronoi_y(2:end) ] ;
            R_j_2  = [mesh(I(i,j)).voronoi_x(end:-1:2) , mesh(I(i,j)).voronoi_y(end:-1:2) , mesh(I(i,j)).voronoi_x(end-1:-1:1) , mesh(I(i,j)).voronoi_y(end-1:-1:1) ]   ;
            JJ     = zeros(length(R_i),1);
            
            for k = 1 : length(R_i(:,1))
                
                for l = 1 : length(R_j_1(:,1))
                    JJ(k) = JJ(k) + sum( ( R_i(k,:)==R_j_1(l,:) )|( R_i(k,:)==R_j_2(l,:) ) );
                end
                
                if (sum(JJ(k))>=6)
                    mesh(i).neighbors(indice) = I(i,j);
                    indice = indice +1;
                end
            end
            
            clear R_j_1 R_j_2 JJ;
        end
        
        
    end


elseif m==3

    x      = C(:,1);
    y      = C(:,2);
    z      = C(:,3);    
    [d]    = give_distance_matrix(x, y,z);
    [~,I]  = sort(d,2,'ascend');
    I      = I(:,2: min(  nb_max_neighbors+1, length(d(:,1)) ) );  
    n_mesh = length(mesh);
    nb_max_neighbors = min(nb_max_neighbors, length(d(:,1)));
    

    for i = 1 : n_mesh

        for j = 1 : length(mesh(i).voronoi_triangular)

                xx = mesh(i).voronoi_triangular(j).x;
                yy = mesh(i).voronoi_triangular(j).y;
                zz = mesh(i).voronoi_triangular(j).z;
                R  = [xx,yy,zz];
                u  = R(1,:)-R(2,:);
                v  = R(1,:)-R(3,:);
                u = u./norm(u);
                v = v./norm(v);
                w  = cross(u,v);
                w  = w./norm(w);
                mesh(i).colinear_face_u(j,:)  =  u;
                mesh(i).colinear_face_v(j,:)  =  v;
                
                mesh(i).normal_face(j,:)  =  w;
                mesh(i).center_faces(j,:) = [ mean(mesh(i).voronoi_triangular(j).x), mean(mesh(i).voronoi_triangular(j).y),...
                                             mean(mesh(i).voronoi_triangular(j).z) ];
        end
    end


    CC = zeros(n_mesh, n_mesh);

    for i = 1 : n_mesh

        for j = 1 :  n_mesh      
      
            w1        = mesh(i).normal_face;
            w2        = mesh(j).normal_face;
            total     = [];
            for kk = 1 : length(w1)
               total  = [total;sum( bsxfun(@times, w1(kk,:), w2),2) ]; 
            end
            total     = abs(total);
            II        = find(total == 1);
            CC(i,j)   = ~isempty(II);
       end
    end

    CC(eye(size(CC))~=0)=0;
    
    II = [1:n_mesh];
    
    for i = 1 : n_mesh
        JJ = find(CC(i,:) == 1);
        mesh(i).neighbors = II(JJ);
    end
    
    

    for  i = 1 : n_mesh

       indice_voisins = mesh(i).neighbors;
       xxi = mesh(i).x;
       yyi = mesh(i).y;
       zzi = mesh(i).z;
       
       for j =1 : length(indice_voisins)
            
           xxj     = mesh(indice_voisins(j)).x;
           yyj     = mesh(indice_voisins(j)).y;
           zzj     = mesh(indice_voisins(j)).z;
           [d]     = give_distance_matrix_2_sets(xxi, yyi,xxj, yyj,zzi,zzj ) ;
           mesh(i).d_min_points_neighbors(j)  = min(d(:)); 
           mesh(i).d_center_neighbors(j) = sqrt( (mesh(i).center_x - mesh(indice_voisins(j)).center_x ).^2 + ...
               (mesh(i).center_y - mesh(indice_voisins(j)).center_y ).^2 +...
               (mesh(i).center_y - mesh(indice_voisins(j)).center_y ).^2);
           mesh(i).ratio_neighbors_close_center(j) = mesh(i).d_min_points_neighbors(j)./mesh(i).d_center_neighbors(j) ;
            
       end
       mesh(i).d_min_points_neighbors_std = std(mesh(i).d_min_points_neighbors);
    end
    
    


end






end