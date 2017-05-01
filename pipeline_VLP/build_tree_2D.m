function [mesh]  = build_tree_2D(tout, parameters, mesh_init, varargin)


format long;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ((nargin ~= 2 ) && (nargin ~= 3 ))
%     
%     if (nargin == 1)
%         disp('too few parameters');
%     else
%         disp('too many parameters');
%     end    
%     
%     mesh = [];
%     
%     return;
%     
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x, y, t,dx , dy] = build_point_map(tout);
    
    if isempty(mesh_init)
        mesh  = init(x,y,t,dx, dy, parameters); 
    else     
        mesh  = mesh_init;
    end
    
    ij = 0;
    
    while(1)
               
        jj = find( [mesh(:).status] == 1);
        
        if isempty(jj)
           break;
        else
            
            indice             = jj(1);
            quadtree_mesh_test = update_quadtree(mesh, indice, parameters);
            status_loc         = [quadtree_mesh_test(indice:indice+3 ).status];
            

            jjj = find( status_loc == 1);
            
            
           
           if ( length(jjj) >= 1)
                critere_reverse = 0;
           else
                critere_reverse = 1; 
           end;
                
      
            
            if critere_reverse ==1
                mesh(indice).status = 0;
            else
                mesh = quadtree_mesh_test;
                clear quadtree_mesh_test;
            end;
            ij =ij + 1;
           
        end;
        
        
    end

    
    mesh                    = postprocessing_position_border(mesh);
    mesh                    = postprocessing_neighbors(mesh);



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh = init(x,y,t,dx, dy,  parameters)

indice = 1;

x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);

mesh(indice).border_left  =  x_min;
mesh(indice).border_right =  x_max;
mesh(indice).border_up    =  y_max;
mesh(indice).border_down  =  y_min;

mesh(indice).center_x    =  (mesh(indice).border_right...
    + mesh(indice).border_left )/2;
mesh(indice).center_y    =  (mesh(indice).border_up ...
    + mesh(indice).border_down )/2;


II = find(  ( x >= mesh(indice).border_left )...
          & ( x <= mesh(indice).border_right)...
          & ( y <= mesh(indice).border_up )...
          & ( y >= mesh(indice).border_down) );

mesh(indice).number_in   =  length(II);

X  = x(II);
Y  = y(II);
T  = t(II);
dX = dx(II);
dY = dy(II);


mesh(indice).x  = X;
mesh(indice).y  = Y;
mesh(indice).t  = T;
mesh(indice).dx = dX;
mesh(indice).dy = dY;
mesh(indice).D         = parameters.D_init;


mesh(indice).status = criterion(mesh, indice, parameters);


    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = criterion(mesh, indice, parameters)

    
    critere_longueur =  generate_critere_longueur(mesh,indice, parameters);
    
    

    if ( ( mesh(indice).number_in >= parameters.number_per_zone) ...
    && ( min(  mesh(indice).border_up - mesh(indice).border_down ...
  , mesh(indice).border_right - mesh(indice).border_left   ) >=  critere_longueur )  )
       
        status = 1;
        
    elseif  mesh(indice).number_in == 0   
        
        status = 2;
                    
    else
        
        status = 0;
        
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh = update_quadtree(mesh, indice, parameters)

n = length(mesh);

for i = 1 : 3;
    mesh( n+i ) = mesh(1);
end;

mesh(indice+4:n+3) =  mesh(indice+1:n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%///////north west//////%
%%%%%%%%%%%%%%%%%%%%%%%%%

%h       = waitbar(0, 'Update Octree ....');
test = mesh(indice);


%waitbar(1/4,h); 


mesh(indice).border_left  =  test.border_left;
mesh(indice).border_right =  (test.border_right + test.border_left)/2;
mesh(indice).border_up    =  test.border_up;
mesh(indice).border_down  =  (test.border_down + test.border_up)/2;


mesh(indice).center_x    =  (mesh(indice).border_right ...
    + mesh(indice).border_left )/2;
mesh(indice).center_y    =  (mesh(indice).border_up    ...
    + mesh(indice).border_down )/2;


II = find(  (test.x >= mesh(indice).border_left  )...
          & (test.x <= mesh(indice).border_right )...
          & (test.y <= mesh(indice).border_up    )...
          & (test.y >= mesh(indice).border_down  ) );

mesh(indice).number_in   =  length(II);


X                   = test.x(II);
Y                   = test.y(II);
T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);

mesh(indice).x      = X;
mesh(indice).y      = Y;
mesh(indice).t      = T;
mesh(indice).dx     = dX;
mesh(indice).dy     = dY;






mesh(indice).status = criterion(mesh, indice, parameters);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% waitbar(2/4,h);
%%%%%%%%%%%%%%%%%%%%%%%%%
%///////north east//////%
%%%%%%%%%%%%%%%%%%%%%%%%%

mesh(indice+1).border_left  =  (test.border_right + test.border_left)/2;
mesh(indice+1).border_right =  test.border_right ;
mesh(indice+1).border_up    =  test.border_up;
mesh(indice+1).border_down  =  (test.border_down + test.border_up)/2;


mesh(indice+1).center_x    =  (mesh(indice+1).border_right ...
    + mesh(indice+1).border_left )/2;
mesh(indice+1).center_y    =  (mesh(indice+1).border_up    ...
    + mesh(indice+1).border_down )/2;


II = find(  (test.x >= mesh(indice+1).border_left  )...
          & (test.x <= mesh(indice+1).border_right )...
          & (test.y <= mesh(indice+1).border_up    )...
          & (test.y >= mesh(indice+1).border_down  ) );

mesh(indice+1).number_in   =  length(II);


X                   = test.x(II);
Y                   = test.y(II);
T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);

mesh(indice+1).x    = X;
mesh(indice+1).y    = Y;
mesh(indice+1).t    = T;
mesh(indice+1).dx   = dX;
mesh(indice+1).dy   = dY;

mesh(indice+1).status = criterion(mesh, indice+1, parameters);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% waitbar(3/4,h);
%%%%%%%%%%%%%%%%%%%%%%%%%
%///////south west//////%
%%%%%%%%%%%%%%%%%%%%%%%%%

mesh(indice+2).border_left  =   test.border_left;
mesh(indice+2).border_right =  (test.border_right + test.border_left)/2; 
mesh(indice+2).border_up    =  (test.border_down  + test.border_up  )/2;
mesh(indice+2).border_down  =   test.border_down;


mesh(indice+2).center_x    =  (mesh(indice+2).border_right ...
    + mesh(indice+2).border_left )/2;
mesh(indice+2).center_y    =  (mesh(indice+2).border_up    ...
    + mesh(indice+2).border_down )/2;


II = find(  (test.x >= mesh(indice+2).border_left  )...
          & (test.x <= mesh(indice+2).border_right )...
          & (test.y <= mesh(indice+2).border_up    )...
          & (test.y >= mesh(indice+2).border_down  ) );

mesh(indice+2).number_in   =  length(II);


X                   = test.x(II);
Y                   = test.y(II);
T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);

mesh(indice+2).x    = X;
mesh(indice+2).y    = Y;
mesh(indice+2).t    = T;
mesh(indice+2).dx   = dX;
mesh(indice+2).dy   = dY;




mesh(indice+2).status = criterion(mesh, indice+2, parameters);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% waitbar(4/4,h);
%%%%%%%%%%%%%%%%%%%%%%%%%
%///////south east//////%
%%%%%%%%%%%%%%%%%%%%%%%%%

mesh(indice+3).border_left  =  (test.border_right + test.border_left)/2;
mesh(indice+3).border_right =  test.border_right ;
mesh(indice+3).border_up    =  (test.border_down  + test.border_up  )/2;
mesh(indice+3).border_down  =  test.border_down;


mesh(indice+3).center_x    =  (mesh(indice+3).border_right ...
    + mesh(indice+3).border_left )/2;
mesh(indice+3).center_y    =  (mesh(indice+3).border_up    ...
    + mesh(indice+3).border_down )/2;


II = find(  (test.x >= mesh(indice+3).border_left  )...
          & (test.x <= mesh(indice+3).border_right )...
          & (test.y <= mesh(indice+3).border_up    )...
          & (test.y >= mesh(indice+3).border_down  ) );

mesh(indice+3).number_in   =  length(II);



X                   = test.x(II);
Y                   = test.y(II);
T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);

mesh(indice+3).x    = X;
mesh(indice+3).y    = Y;
mesh(indice+3).t    = T;
mesh(indice+3).dx   = dX;
mesh(indice+3).dy   = dY;




mesh(indice+3).status = criterion(mesh, indice+3, parameters);

% close(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh = postprocessing_position_border(mesh)

n = length(mesh);

% h       = waitbar(0, 'Check the borders ....');

for i =1 : n
%    waitbar(i/n,h);
    
   mesh(i).box(1).x     = [mesh(i).border_left   mesh(i).border_left];
   mesh(i).box(1).y     = [mesh(i).border_down   mesh(i).border_up];
   
   mesh(i).box(2).x     = [mesh(i).border_right  mesh(i).border_right];
   mesh(i).box(2).y     = [mesh(i).border_down   mesh(i).border_up];
   
   mesh(i).box(3).x     = [mesh(i).border_left   mesh(i).border_right];
   mesh(i).box(3).y     = [mesh(i).border_down   mesh(i).border_down];
   
   mesh(i).box(4).x     = [mesh(i).border_left   mesh(i).border_right];
   mesh(i).box(4).y     = [mesh(i).border_up     mesh(i).border_up];
   
   mesh(i).cornerx(1)   = mesh(i).border_left     ;
   mesh(i).cornery(1)   = mesh(i).border_up       ;
   
   mesh(i).cornerx(2)   = mesh(i).border_left     ;
   mesh(i).cornery(2)   = mesh(i).border_down     ;
   
   mesh(i).cornerx(3)   = mesh(i).border_right    ;
   mesh(i).cornery(3)   = mesh(i).border_up       ;
   
   mesh(i).cornerx(4)   = mesh(i).border_right    ;
   mesh(i).cornery(4)   = mesh(i).border_down     ;
   
   mesh(i).tree_x       = [mesh(i).cornerx(1), mesh(i).cornerx(2), mesh(i).cornerx(4),  mesh(i).cornerx(3) , mesh(i).cornerx(1)];
   mesh(i).tree_y       = [mesh(i).cornery(1), mesh(i).cornery(2), mesh(i).cornery(4),  mesh(i).cornery(3) , mesh(i).cornery(1)];
   
   
   
   
end

% close(h);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh = postprocessing_neighbors(mesh)


n = length(mesh);
for i =1 : n
    mesh(i).numero = i;
end;

corners = zeros(n,4,2);


% h       = waitbar(0, 'Find Corners ....');


for i = 1 : n;
%     waitbar(i/n,h);
    corners(i,1,1) = mesh(i).border_left;
    corners(i,1,2) = mesh(i).border_up;

    corners(i,2,1) = mesh(i).border_left;
    corners(i,2,2) = mesh(i).border_down;

    corners(i,3,1) = mesh(i).border_right;
    corners(i,3,2) = mesh(i).border_up;

    corners(i,4,1) = mesh(i).border_right;
    corners(i,4,2) = mesh(i).border_down;

end;
% close(h);

JJ(1,:) = [2,3,4];
JJ(2,:) = [1,3,4];
JJ(3,:) = [1,2,4];
JJ(4,:) = [1,2,3];

% h       = waitbar(0, 'Search for Neighbors ....');

for i = 1 : n;
%     waitbar(i/n,h);
    II = [];
    for k = 1 : n;
        for j = 1 : 4;
            for l = 1:3
                if ( ( corners(i,j,1) == corners(k,JJ(j,l),1)  ) && ( corners(i,j,2) == corners(k,JJ(j,l),2)  )  )
                    II =[II; k];
                end;
            end;
        end;
    end;

    

    for k = 1 : n;

            if (  ( ( mesh(i).box(1).y(1)  >=  mesh(k).box(2).y(1)  ) && (   mesh(i).box(1).y(2)  <=  mesh(k).box(2).y(2)    ) && ( mesh(i).box(1).x(1) == mesh(k).box(2).x(1) ) ) ...
            ||    ( ( mesh(i).box(1).y(1)  <=  mesh(k).box(2).y(1)  ) && (   mesh(i).box(1).y(2)  >=  mesh(k).box(2).y(2)    ) && ( mesh(i).box(1).x(1) == mesh(k).box(2).x(1) ) ) )
            II =[II; k];
            end
       
             if (  ( ( mesh(i).box(2).y(1)  >=  mesh(k).box(1).y(1)  ) && (   mesh(i).box(2).y(2)  <=  mesh(k).box(1).y(2)    ) ) && (mesh(i).box(2).x(1)  ==  mesh(k).box(1).x(1) ) ...
            ||     ( ( mesh(i).box(2).y(1)  <=  mesh(k).box(1).y(1)  ) && (   mesh(i).box(2).y(2)  >=  mesh(k).box(1).y(2)    ) ) && (mesh(i).box(2).x(1)  ==  mesh(k).box(1).x(1) ))
            II =[II; k];
             end
            
            if (  ( ( mesh(i).box(3).x(1)  >=  mesh(k).box(4).x(1)  ) && (   mesh(i).box(3).x(2)  <=  mesh(k).box(4).x(2)    ) ) && (mesh(i).box(3).y(1)  ==  mesh(k).box(4).y(1))...
            ||    ( ( mesh(i).box(3).x(1)  <=  mesh(k).box(4).x(1)  ) && (   mesh(i).box(3).x(2)  >=  mesh(k).box(4).x(2)    ) ) && (mesh(i).box(3).y(1)  ==  mesh(k).box(4).y(1)) )
            II =[II; k];
             end
             
             if (  ( ( mesh(i).box(4).x(1)  >=  mesh(k).box(3).x(1)  ) && (   mesh(i).box(4).x(2)  <=  mesh(k).box(3).x(2)    ) ) &&  (mesh(i).box(4).y(1)  ==  mesh(k).box(3).y(1) )...
            ||     ( ( mesh(i).box(4).x(1)  <=  mesh(k).box(3).x(1)  ) && (   mesh(i).box(4).x(2)  >=  mesh(k).box(3).x(2)    ) ) &&  (mesh(i).box(4).y(1)  ==  mesh(k).box(3).y(1) ) )
            II =[II; k];
            end
        
    end;
    II              = unique(II);
    mesh(i).neighbors = II;
end;
% close(h);



% h       = waitbar(0, 'Organize neighbors ....');
for i = 1 : n;
%     waitbar(i/n,h);
    kk = mesh(i).neighbors;

    mesh(i).plus_x          = [];
    mesh(i).plus_x_index    = [];
    mesh(i).minus_x         = [];
    mesh(i).minus_x_index   = [];
    mesh(i).plus_y          = [];
    mesh(i).plus_y_index    = [];
    mesh(i).minus_y         = [];
    mesh(i).minus_y_index   = [];


    %% 1
    xxi = mesh(i).box(1).x;
    yyi = mesh(i).box(1).y;

    for j = 1 : length(kk);
        kkkk = kk(j);
        
        xxj = mesh(kkkk).box(2).x;
        yyj = mesh(kkkk).box(2).y;
        
        II  = inpolygon( xxi, yyi, xxj, yyj);
        JJ  = inpolygon( xxj, yyj, xxi, yyi);
        
            if ((sum(II)==2)||(sum(JJ)==2))
                mesh(i).minus_x = [mesh(i).minus_x; mesh(kkkk).center_x];
                mesh(i).minus_x_index = [mesh(i).minus_x_index; kkkk];
            end
        
    end

    %% 2
    xxi = mesh(i).box(2).x;
    yyi = mesh(i).box(2).y;

    for j = 1 : length(kk);
        kkkk = kk(j);
        
        xxj = mesh(kkkk).box(1).x;
        yyj = mesh(kkkk).box(1).y;
        
        II  = inpolygon( xxi, yyi, xxj, yyj);
        JJ  = inpolygon( xxj, yyj, xxi, yyi);
        
            if ((sum(II)==2)||(sum(JJ)==2))
                mesh(i).plus_x = [mesh(i).plus_x; mesh(kkkk).center_x];
                mesh(i).plus_x_index = [mesh(i).plus_x_index; kkkk];
            end
        
    end
    
    

    %% 3
    xxi = mesh(i).box(3).x;
    yyi = mesh(i).box(3).y;

    for j = 1 : length(kk);
        kkkk = kk(j);
        
        xxj = mesh(kkkk).box(4).x;
        yyj = mesh(kkkk).box(4).y;
        
        II  = inpolygon( xxi, yyi, xxj, yyj);
        JJ  = inpolygon( xxj, yyj, xxi, yyi);
        
            if ((sum(II)==2)||(sum(JJ)==2))
                mesh(i).minus_y = [mesh(i).minus_y; mesh(kkkk).center_y];
                mesh(i).minus_y_index = [mesh(i).minus_y_index; kkkk];
            end
        
    end
    

    %% 4
    xxi = mesh(i).box(4).x;
    yyi = mesh(i).box(4).y;

    for j = 1 : length(kk);
        kkkk = kk(j);
        
        xxj = mesh(kkkk).box(3).x;
        yyj = mesh(kkkk).box(3).y;
        
        II  = inpolygon( xxi, yyi, xxj, yyj);
        JJ  = inpolygon( xxj, yyj, xxi, yyi);
        
            if ((sum(II)==2)||(sum(JJ)==2))
                mesh(i).plus_y = [mesh(i).plus_y; mesh(kkkk).center_y];
                mesh(i).plus_y_index = [mesh(i).plus_y_index; kkkk];
            end
        
    end



    
end
% close(h);
    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mesh, movie_per_frame] = postprocessing_trajectories(mesh, x, y, t, movie_per_frame )


for i = 1 : length(mesh)
     
    II      = find(       ( x >= mesh(i).border_left  )...
                        & ( x <= mesh(i).border_right )...
                        & ( y <= mesh(i).border_up    )...
                        & ( y >= mesh(i).border_down  ) );

    X                                   = x(II);
    Y                                   = y(II);
    T                                   = t(II);
    TT                                  = unique(T);
    
    for jj = 1 : length(TT)
        tt_loc                          = TT(jj);
        jjj                             = find( (T == tt_loc) );
        xx                              = X(jjj);
        yy                              = Y(jjj);
        mesh(i).movie_per_frame(jj).t   = tt_loc;
        mesh(i).movie_per_frame(jj).x   = xx;
        mesh(i).movie_per_frame(jj).y   = yy;
    end;
    
end;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = generate_critere_longueur(mesh,indice, parameters)

res                 = 1.*sqrt(2*mesh(indice).D*parameters.dt_theo);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

