function [mesh]  = build_tree_3D(tout, parameters, mesh_init, varargin)



format long;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [x, y, t,dx , dy] = build_point_map(tout);
    if isempty(mesh_init) == 2
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
            octree_mesh_test = update_octree(mesh, indice, parameters);
           % critere_reverse    = 0;
            status_loc         = [octree_mesh_test(indice:indice+3 ).status];
            
            %iii = find( status_loc == 0);
            jjj = find( status_loc == 1);
            
            
           
           if ( length(jjj) >= 1)
                critere_reverse = 0;
           else
                critere_reverse = 1; 
           end;
                
      
            
            if critere_reverse ==1
                mesh(indice).status = 0;
            else
                mesh = octree_mesh_test;
                clear octree_mesh_test;
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
function mesh = init(x,y,z,t,dx, dy,dz,  parameters)

indice = 1;

x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);
z_min = min(z);
z_max = max(z);



mesh(indice).border_left     =  x_min;
mesh(indice).border_right    =  x_max;
mesh(indice).border_up       =  y_max;
mesh(indice).border_down     =  y_min;
mesh(indice).border_north    =  z_max;
mesh(indice).border_south    =  z_min;


mesh(indice).center_x    =  (mesh(indice).border_right...
    + mesh(indice).border_left )/2;
mesh(indice).center_y    =  (mesh(indice).border_up ...
    + mesh(indice).border_down )/2;
mesh(indice).center_z    =  (mesh(indice).border_north ...
    + mesh(indice).border_south )/2;



II = find(  ( x >= mesh(indice).border_left )...
          & ( x <= mesh(indice).border_right)...
          & ( y <= mesh(indice).border_up )...
          & ( y >= mesh(indice).border_down) ...
          & ( z <= mesh(indice).border_north )...
          & ( z >= mesh(indice).border_south) ...
          );

mesh(indice).number_in   =  length(II);

X  = x(II);
Y  = y(II);
Z  = z(II);
T  = t(II);
dX = dx(II);
dY = dy(II);
dZ = dz(II);



mesh(indice).x  = X;
mesh(indice).y  = Y;
mesh(indice).z  = Z;

mesh(indice).t  = T;
mesh(indice).dx = dX;
mesh(indice).dy = dY;
mesh(indice).dz = dZ;

mesh(indice).D         = parameters.D_init;


mesh(indice).status = criterion(mesh, indice, parameters);


    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = criterion(mesh, indice, parameters)

    
    critere_longueur =  generate_critere_longueur(mesh,indice, parameters);
    
    

    if ( ( mesh(indice).number_in >= parameters.number_per_zone) ...
    && ( min( min(  mesh(indice).border_up - mesh(indice).border_down ...
  , mesh(indice).border_right - mesh(indice).border_left) ,...
  mesh(indice).border_north - mesh(indice).border_south  ) >=  critere_longueur )  )
       
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
function mesh = update_octree(mesh, indice, parameters)

n = length(mesh);
% 
% for i = 1 : 3;
%     mesh( n+i ) = mesh(1);
% end;
% 
% mesh(indice+4:n+3) =  mesh(indice+1:n);

for i = 1 : 7;
    mesh( n+i ) = mesh(1);
end;

mesh(indice+8:n+7) =  mesh(indice+1:n);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h       = waitbar(0, 'Update Octree ....');
test = mesh(indice);


%waitbar(1/8,h);


%%%%%%%%%%%%%%%%%%%%%%%%%
%///////up north west//////%
%%%%%%%%%%%%%%%%%%%%%%%%%






mesh(indice).border_left     =  test.border_left;
mesh(indice).border_right    =  (test.border_right + test.border_left)/2;
mesh(indice).border_up       =  test.border_up;
mesh(indice).border_down     =  (test.border_down + test.border_up)/2;
mesh(indice).border_north    =  test.border_north;
mesh(indice).border_south    =  (test.border_north + test.border_south)/2;



mesh(indice).center_x    =  (mesh(indice).border_right ...
    + mesh(indice).border_left )/2;
mesh(indice).center_y    =  (mesh(indice).border_up    ...
    + mesh(indice).border_down )/2;
mesh(indice).center_z    =  (mesh(indice).border_north    ...
    + mesh(indice).border_south )/2;




      II = find(  ( test.x >= mesh(indice).border_left )...
          & ( test.x <= mesh(indice).border_right)...
          & ( test.y <= mesh(indice).border_up )...
          & ( test.y >= mesh(indice).border_down) ...
          & ( test.z <= mesh(indice).border_north )...
          & ( test.z >= mesh(indice).border_south) ...
          );
      
      
mesh(indice).number_in   =  length(II);


X                   = test.x(II);
Y                   = test.y(II);
Z                   = test.z(II);
T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);
dZ                  = test.dz(II);


mesh(indice).x      = X;
mesh(indice).y      = Y;
mesh(indice).z      = Z;

mesh(indice).t      = T;
mesh(indice).dx     = dX;
mesh(indice).dy     = dY;
mesh(indice).dz     = dZ;



mesh(indice).status = criterion(mesh, indice, parameters);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% waitbar(2/8,h);
%%%%%%%%%%%%%%%%%%%%%%%%%
%///////up north east//////%
%%%%%%%%%%%%%%%%%%%%%%%%%

mesh(indice+1).border_left  =  (test.border_right + test.border_left)/2;
mesh(indice+1).border_right =  test.border_right ;
mesh(indice+1).border_up    =  test.border_up;
mesh(indice+1).border_down  =  (test.border_down + test.border_up)/2;
mesh(indice+1).border_north    =  test.border_north;
mesh(indice+1).border_south    =  (test.border_north + test.border_south)/2;

mesh(indice+1).center_x    =  (mesh(indice+1).border_right ...
    + mesh(indice+1).border_left )/2;
mesh(indice+1).center_y    =  (mesh(indice+1).border_up    ...
    + mesh(indice+1).border_down )/2;
mesh(indice+1).center_z    =  (mesh(indice+1).border_north    ...
    + mesh(indice+1).border_south )/2;

      
  II = find(  ( test.x >= mesh(indice+1).border_left )...
          & ( test.x <= mesh(indice+1).border_right)...
          & ( test.y <= mesh(indice+1).border_up )...
          & ( test.y >= mesh(indice+1).border_down) ...
          & ( test.z <= mesh(indice+1).border_north )...
          & ( test.z >= mesh(indice+1).border_south) ...
          );      
      

mesh(indice+1).number_in   =  length(II);


X                   = test.x(II);
Y                   = test.y(II);
Z                   = test.z(II);
T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);
dZ                  = test.dz(II);


mesh(indice+1).x    = X;
mesh(indice+1).y    = Y;
mesh(indice+1).z    = Z;
mesh(indice+1).t    = T;
mesh(indice+1).dx   = dX;
mesh(indice+1).dy   = dY;
mesh(indice+1).dz   = dZ;


mesh(indice+1).status = criterion(mesh, indice+1, parameters);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% waitbar(3/8,h);
%%%%%%%%%%%%%%%%%%%%%%%%%
%///////up south west//////%
%%%%%%%%%%%%%%%%%%%%%%%%%

mesh(indice+2).border_left  =   test.border_left;
mesh(indice+2).border_right =  (test.border_right + test.border_left)/2; 
mesh(indice+2).border_up    =  (test.border_down  + test.border_up  )/2;
mesh(indice+2).border_down  =   test.border_down;
mesh(indice+2).border_north    =  test.border_north;
mesh(indice+2).border_south    =  (test.border_north + test.border_south)/2;

mesh(indice+2).center_x    =  (mesh(indice+2).border_right ...
    + mesh(indice+2).border_left )/2;
mesh(indice+2).center_y    =  (mesh(indice+2).border_up    ...
    + mesh(indice+2).border_down )/2;
mesh(indice+2).center_z    =  (mesh(indice+2).border_north    ...
    + mesh(indice+2).border_south )/2;

  II = find(  ( test.x >= mesh(indice+2).border_left )...
          & ( test.x <= mesh(indice+2).border_right)...
          & ( test.y <= mesh(indice+2).border_up )...
          & ( test.y >= mesh(indice+2).border_down) ...
          & ( test.z <= mesh(indice+2).border_north )...
          & ( test.z >= mesh(indice+2).border_south) ...
          );        
      
      
mesh(indice+2).number_in   =  length(II);


X                   = test.x(II);
Y                   = test.y(II);
Z                   = test.z(II);
T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);
dZ                  = test.dz(II);


mesh(indice+2).x    = X;
mesh(indice+2).y    = Y;
mesh(indice+2).z    = Z;

mesh(indice+2).t    = T;
mesh(indice+2).dx   = dX;
mesh(indice+2).dy   = dY;
mesh(indice+2).dz   = dZ;





mesh(indice+2).status = criterion(mesh, indice+2, parameters);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% waitbar(4/8,h);
%%%%%%%%%%%%%%%%%%%%%%%%%
%///////up south east//////%
%%%%%%%%%%%%%%%%%%%%%%%%%

mesh(indice+3).border_left  =  (test.border_right + test.border_left)/2;
mesh(indice+3).border_right =  test.border_right ;
mesh(indice+3).border_up    =  (test.border_down  + test.border_up  )/2;
mesh(indice+3).border_down  =  test.border_down;
mesh(indice+3).border_north    =  test.border_north;
mesh(indice+3).border_south    =  (test.border_north + test.border_south)/2;

mesh(indice+3).center_x    =  (mesh(indice+3).border_right ...
    + mesh(indice+3).border_left )/2;
mesh(indice+3).center_y    =  (mesh(indice+3).border_up    ...
    + mesh(indice+3).border_down )/2;
mesh(indice+3).center_z    =  (mesh(indice+3).border_north    ...
    + mesh(indice+3).border_south )/2;
      
      
  II = find(  ( test.x >= mesh(indice+3).border_left )...
          & ( test.x <= mesh(indice+3).border_right)...
          & ( test.y <= mesh(indice+3).border_up )...
          & ( test.y >= mesh(indice+3).border_down) ...
          & ( test.z <= mesh(indice+3).border_north )...
          & ( test.z >= mesh(indice+3).border_south) ...
          );   

mesh(indice+3).number_in   =  length(II);



X                   = test.x(II);
Y                   = test.y(II);
Z                   = test.z(II);

T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);
dZ                  = test.dz(II);


mesh(indice+3).x    = X;
mesh(indice+3).y    = Y;
mesh(indice+3).z    = Z;

mesh(indice+3).t    = T;
mesh(indice+3).dx   = dX;
mesh(indice+3).dy   = dY;
mesh(indice+3).dz   = dZ;





mesh(indice+3).status = criterion(mesh, indice+3, parameters);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














% waitbar(5/8,h);
%%%%%%%%%%%%%%%%%%%%%%%%%
%///////down north west//////%
%%%%%%%%%%%%%%%%%%%%%%%%%

mesh(indice+4).border_left     =  test.border_left;
mesh(indice+4).border_right    =  (test.border_right + test.border_left)/2;
mesh(indice+4).border_up       =  test.border_up;
mesh(indice+4).border_down     =  (test.border_down + test.border_up)/2;
mesh(indice+4).border_north    =  (test.border_north + test.border_south)/2;
mesh(indice+4).border_south    =  ( test.border_south);



mesh(indice+4).center_x    =  (mesh(indice+4).border_right ...
    + mesh(indice+4).border_left )/2;
mesh(indice+4).center_y    =  (mesh(indice+4).border_up    ...
    + mesh(indice+4).border_down )/2;
mesh(indice+4).center_z    =  (mesh(indice+4).border_north    ...
    + mesh(indice+4).border_south )/2;




      II = find(  ( test.x >= mesh(indice+4).border_left )...
          & ( test.x <= mesh(indice+4).border_right)...
          & ( test.y <= mesh(indice+4).border_up )...
          & ( test.y >= mesh(indice+4).border_down) ...
          & ( test.z <= mesh(indice+4).border_north )...
          & ( test.z >= mesh(indice+4).border_south) ...
          );
      
      
mesh(indice+4).number_in   =  length(II);


X                   = test.x(II);
Y                   = test.y(II);
Z                   = test.z(II);
T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);
dZ                  = test.dz(II);


mesh(indice+4).x      = X;
mesh(indice+4).y      = Y;
mesh(indice+4).z      = Z;

mesh(indice+4).t      = T;
mesh(indice+4).dx     = dX;
mesh(indice+4).dy     = dY;
mesh(indice+4).dz     = dZ;



mesh(indice+4).status = criterion(mesh, indice+4, parameters);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% waitbar(6/8,h);
%%%%%%%%%%%%%%%%%%%%%%%%%
%///////down north east//////%
%%%%%%%%%%%%%%%%%%%%%%%%%

mesh(indice+5).border_left  =  (test.border_right + test.border_left)/2;
mesh(indice+5).border_right =  test.border_right ;
mesh(indice+5).border_up    =  test.border_up;
mesh(indice+5).border_down  =  (test.border_down + test.border_up)/2;
mesh(indice+5).border_north    =  (test.border_north + test.border_south)/2;
mesh(indice+5).border_south    =  ( test.border_south);

mesh(indice+5).center_x    =  (mesh(indice+5).border_right ...
    + mesh(indice+5).border_left )/2;
mesh(indice+5).center_y    =  (mesh(indice+5).border_up    ...
    + mesh(indice+5).border_down )/2;
mesh(indice+5).center_z    =  (mesh(indice+5).border_north    ...
    + mesh(indice+5).border_south )/2;

      
  II = find(  ( test.x >= mesh(indice+5).border_left )...
          & ( test.x <= mesh(indice+5).border_right)...
          & ( test.y <= mesh(indice+5).border_up )...
          & ( test.y >= mesh(indice+5).border_down) ...
          & ( test.z <= mesh(indice+5).border_north )...
          & ( test.z >= mesh(indice+5).border_south) ...
          );      
      

mesh(indice+5).number_in   =  length(II);


X                   = test.x(II);
Y                   = test.y(II);
Z                   = test.z(II);
T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);
dZ                  = test.dz(II);


mesh(indice+5).x    = X;
mesh(indice+5).y    = Y;
mesh(indice+5).z    = Z;
mesh(indice+5).t    = T;
mesh(indice+5).dx   = dX;
mesh(indice+5).dy   = dY;
mesh(indice+5).dz   = dZ;


mesh(indice+5).status = criterion(mesh, indice+5, parameters);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% waitbar(7/8,h);
%%%%%%%%%%%%%%%%%%%%%%%%%
%///////down south west//////%
%%%%%%%%%%%%%%%%%%%%%%%%%

mesh(indice+6).border_left  =   test.border_left;
mesh(indice+6).border_right =  (test.border_right + test.border_left)/2; 
mesh(indice+6).border_up    =  (test.border_down  + test.border_up  )/2;
mesh(indice+6).border_down  =   test.border_down;
mesh(indice+6).border_north    =  (test.border_north + test.border_south)/2;
mesh(indice+6).border_south    =  ( test.border_south);

mesh(indice+6).center_x    =  (mesh(indice+6).border_right ...
    + mesh(indice+6).border_left )/2;
mesh(indice+6).center_y    =  (mesh(indice+6).border_up    ...
    + mesh(indice+6).border_down )/2;
mesh(indice+6).center_z    =  (mesh(indice+6).border_north    ...
    + mesh(indice+6).border_south )/2;

  II = find(  ( test.x >= mesh(indice+6).border_left )...
          & ( test.x <= mesh(indice+6).border_right)...
          & ( test.y <= mesh(indice+6).border_up )...
          & ( test.y >= mesh(indice+6).border_down) ...
          & ( test.z <= mesh(indice+6).border_north )...
          & ( test.z >= mesh(indice+6).border_south) ...
          );        
      
      
mesh(indice+6).number_in   =  length(II);


X                   = test.x(II);
Y                   = test.y(II);
Z                   = test.z(II);
T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);
dZ                  = test.dz(II);


mesh(indice+6).x    = X;
mesh(indice+6).y    = Y;
mesh(indice+6).z    = Z;
mesh(indice+6).t    = T;
mesh(indice+6).dx   = dX;
mesh(indice+6).dy   = dY;
mesh(indice+6).dz   = dZ;





mesh(indice+6).status = criterion(mesh, indice+6, parameters);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% waitbar(8/8,h);
%%%%%%%%%%%%%%%%%%%%%%%%%
%///////down south east//////%
%%%%%%%%%%%%%%%%%%%%%%%%%

mesh(indice+7).border_left  =  (test.border_right + test.border_left)/2;
mesh(indice+7).border_right =  test.border_right ;
mesh(indice+7).border_up    =  (test.border_down  + test.border_up  )/2;
mesh(indice+7).border_down  =  test.border_down;
mesh(indice+7).border_north    =  (test.border_north + test.border_south)/2;
mesh(indice+7).border_south    =  ( test.border_south);

mesh(indice+7).center_x    =  (mesh(indice+7).border_right ...
    + mesh(indice+7).border_left )/2;
mesh(indice+7).center_y    =  (mesh(indice+7).border_up    ...
    + mesh(indice+7).border_down )/2;
mesh(indice+7).center_z    =  (mesh(indice+7).border_north    ...
    + mesh(indice+7).border_south )/2;
      
      
  II = find(  ( test.x >= mesh(indice+7).border_left )...
          & ( test.x <= mesh(indice+7).border_right)...
          & ( test.y <= mesh(indice+7).border_up )...
          & ( test.y >= mesh(indice+7).border_down) ...
          & ( test.z <= mesh(indice+7).border_north )...
          & ( test.z >= mesh(indice+7).border_south) ...
          );   

mesh(indice+7).number_in   =  length(II);



X                   = test.x(II);
Y                   = test.y(II);
Z                   = test.z(II);

T                   = test.t(II);
dX                  = test.dx(II);
dY                  = test.dy(II);
dZ                  = test.dz(II);


mesh(indice+7).x    = X;
mesh(indice+7).y    = Y;
mesh(indice+7).z    = Z;

mesh(indice+7).t    = T;
mesh(indice+7).dx   = dX;
mesh(indice+7).dy   = dY;
mesh(indice+7).dz   = dZ;





mesh(indice+7).status = criterion(mesh, indice+7, parameters);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% close(h);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh = postprocessing_position_border(mesh)

n = length(mesh);

% h       = waitbar(0, 'Check the borders ....');

for i =1 : n
%    waitbar(i/n,h);
    %%
   mesh(i).box(1).x     = [mesh(i).border_left    mesh(i).border_right];
   mesh(i).box(1).y     = [mesh(i).border_down    mesh(i).border_down];
   mesh(i).box(1).z     = [mesh(i).border_north   mesh(i).border_north];
   
   mesh(i).box(2).x     = [mesh(i).border_left    mesh(i).border_right];
   mesh(i).box(2).y     = [mesh(i).border_up      mesh(i).border_up];
   mesh(i).box(2).z     = [mesh(i).border_north   mesh(i).border_north];
   
   mesh(i).box(3).x     = [mesh(i).border_left    mesh(i).border_right];
   mesh(i).box(3).y     = [mesh(i).border_down    mesh(i).border_down];
   mesh(i).box(3).z     = [mesh(i).border_south   mesh(i).border_south];
   
   mesh(i).box(4).x     = [mesh(i).border_left    mesh(i).border_right];
   mesh(i).box(4).y     = [mesh(i).border_up      mesh(i).border_up];
   mesh(i).box(4).z     = [mesh(i).border_south   mesh(i).border_south];
   
   
   %%
   mesh(i).box(5).x     = [mesh(i).border_left    mesh(i).border_left];
   mesh(i).box(5).y     = [mesh(i).border_down    mesh(i).border_up];
   mesh(i).box(5).z     = [mesh(i).border_north   mesh(i).border_north];
   
   mesh(i).box(6).x     = [mesh(i).border_right    mesh(i).border_right];
   mesh(i).box(6).y     = [mesh(i).border_down    mesh(i).border_up];
   mesh(i).box(6).z     = [mesh(i).border_north   mesh(i).border_north];
   
   mesh(i).box(7).x     = [mesh(i).border_left    mesh(i).border_left];
   mesh(i).box(7).y     = [mesh(i).border_down    mesh(i).border_up];
   mesh(i).box(7).z     = [mesh(i).border_south   mesh(i).border_south];
   
   mesh(i).box(8).x     = [mesh(i).border_right    mesh(i).border_right];
   mesh(i).box(8).y     = [mesh(i).border_down    mesh(i).border_up];
   mesh(i).box(8).z     = [mesh(i).border_south   mesh(i).border_south];
   
   %%
   mesh(i).box(9).x     = [mesh(i).border_left    mesh(i).border_left];
   mesh(i).box(9).y     = [mesh(i).border_down    mesh(i).border_down];
   mesh(i).box(9).z     = [mesh(i).border_south   mesh(i).border_north];
   
   mesh(i).box(10).x     = [mesh(i).border_left    mesh(i).border_left];
   mesh(i).box(10).y     = [mesh(i).border_up    mesh(i).border_up];
   mesh(i).box(10).z     = [mesh(i).border_south   mesh(i).border_north];
   
   mesh(i).box(11).x     = [mesh(i).border_right    mesh(i).border_right];
   mesh(i).box(11).y     = [mesh(i).border_down    mesh(i).border_down];
   mesh(i).box(11).z     = [mesh(i).border_south   mesh(i).border_north];
   
   mesh(i).box(12).x     = [mesh(i).border_right    mesh(i).border_right];
   mesh(i).box(12).y     = [mesh(i).border_up    mesh(i).border_up];
   mesh(i).box(12).z     = [mesh(i).border_south   mesh(i).border_north];
   
   
   %%
   
   
   mesh(i).cornerx(1)   = mesh(i).border_left     ;
   mesh(i).cornery(1)   = mesh(i).border_down       ;
   mesh(i).cornerz(1)   = mesh(i).border_south       ;
   
   mesh(i).cornerx(2)   = mesh(i).border_left     ;
   mesh(i).cornery(2)   = mesh(i).border_down       ;
   mesh(i).cornerz(2)   = mesh(i).border_north       ;
   
   mesh(i).cornerx(3)   = mesh(i).border_left     ;
   mesh(i).cornery(3)   = mesh(i).border_up       ;
   mesh(i).cornerz(3)   = mesh(i).border_south       ;
   
   mesh(i).cornerx(4)   = mesh(i).border_left     ;
   mesh(i).cornery(4)   = mesh(i).border_up       ;
   mesh(i).cornerz(4)   = mesh(i).border_north       ;
   
   mesh(i).cornerx(5)   = mesh(i).border_right     ;
   mesh(i).cornery(5)   = mesh(i).border_down       ;
   mesh(i).cornerz(5)   = mesh(i).border_south       ;
   
   mesh(i).cornerx(6)   = mesh(i).border_right     ;
   mesh(i).cornery(6)   = mesh(i).border_down       ;
   mesh(i).cornerz(6)   = mesh(i).border_north       ;
   
   mesh(i).cornerx(7)   = mesh(i).border_right    ;
   mesh(i).cornery(7)   = mesh(i).border_up       ;
   mesh(i).cornerz(7)   = mesh(i).border_south       ;
   
   mesh(i).cornerx(8)   = mesh(i).border_right     ;
   mesh(i).cornery(8)   = mesh(i).border_up       ;
   mesh(i).cornerz(8)   = mesh(i).border_north       ;
   
%    mesh(i).tree_x       = [mesh(i).cornerx(1), mesh(i).cornerx(2), mesh(i).cornerx(4),  mesh(i).cornerx(3) , mesh(i).cornerx(1)];
%    mesh(i).tree_y       = [mesh(i).cornery(1), mesh(i).cornery(2), mesh(i).cornery(4),  mesh(i).cornery(3) , mesh(i).cornery(1)];
%    mesh(i).tree_z       = [mesh(i).cornery(1), mesh(i).cornery(2), mesh(i).cornery(4),  mesh(i).cornery(3) , mesh(i).cornery(1)];
%    
%    
   
   %% faces
   
   %% 1
   mesh(i).faces(1).xmin = mesh(i).border_left;
   mesh(i).faces(1).xmax = mesh(i).border_right;
   mesh(i).faces(1).ymin = mesh(i).border_down;
   mesh(i).faces(1).ymax = mesh(i).border_up ;
   mesh(i).faces(1).zmin = mesh(i).border_north ;
   mesh(i).faces(1).zmax = mesh(i).border_north ;
   mesh(i).faces(1).square_x = [mesh(i).faces(1).xmin ;mesh(i).faces(1).xmax; mesh(i).faces(1).xmax; mesh(i).faces(1).xmin; mesh(i).faces(1).xmin ;];
   mesh(i).faces(1).square_y = [mesh(i).faces(1).ymin ;mesh(i).faces(1).ymin; mesh(i).faces(1).ymax; mesh(i).faces(1).ymax; mesh(i).faces(1).ymin; ];
   mesh(i).faces(1).square_z = [mesh(i).faces(1).zmin ;mesh(i).faces(1).zmin ;mesh(i).faces(1).zmin ;mesh(i).faces(1).zmin ;mesh(i).faces(1).zmin ;];
   
   
   
   %% 2
   mesh(i).faces(2).xmin = mesh(i).border_left;
   mesh(i).faces(2).xmax = mesh(i).border_right;
   mesh(i).faces(2).ymin = mesh(i).border_down;
   mesh(i).faces(2).ymax = mesh(i).border_up ;
   mesh(i).faces(2).zmin = mesh(i).border_south ;
   mesh(i).faces(2).zmax = mesh(i).border_south ;
   mesh(i).faces(2).square_x = [mesh(i).faces(2).xmin ;mesh(i).faces(2).xmax; mesh(i).faces(2).xmax; mesh(i).faces(2).xmin; mesh(i).faces(2).xmin ;];
   mesh(i).faces(2).square_y = [mesh(i).faces(2).ymin ;mesh(i).faces(2).ymin; mesh(i).faces(2).ymax; mesh(i).faces(2).ymax; mesh(i).faces(2).ymin; ];
   mesh(i).faces(2).square_z = [mesh(i).faces(2).zmin ;mesh(i).faces(2).zmin ;mesh(i).faces(2).zmin ;mesh(i).faces(2).zmin ;mesh(i).faces(2).zmin ;];
   
   
    %% 3
   mesh(i).faces(3).xmin = mesh(i).border_right;
   mesh(i).faces(3).xmax = mesh(i).border_right;
   mesh(i).faces(3).ymin = mesh(i).border_down;
   mesh(i).faces(3).ymax = mesh(i).border_up ;
   mesh(i).faces(3).zmin = mesh(i).border_south ;
   mesh(i).faces(3).zmax = mesh(i).border_north ;
   mesh(i).faces(3).square_y = [mesh(i).faces(3).ymin ;mesh(i).faces(3).ymax; mesh(i).faces(3).ymax; mesh(i).faces(3).ymin; mesh(i).faces(3).ymin ;];
   mesh(i).faces(3).square_z = [mesh(i).faces(3).zmin ;mesh(i).faces(3).zmin; mesh(i).faces(3).zmax; mesh(i).faces(3).zmax; mesh(i).faces(3).zmin; ];
   mesh(i).faces(3).square_x = [mesh(i).faces(3).xmin ;mesh(i).faces(3).xmax; mesh(i).faces(3).xmax; mesh(i).faces(3).xmin; mesh(i).faces(3).xmin ;];
   
   
    %% 4
   mesh(i).faces(4).xmin = mesh(i).border_left;
   mesh(i).faces(4).xmax = mesh(i).border_right;
   mesh(i).faces(4).ymin = mesh(i).border_up;
   mesh(i).faces(4).ymax = mesh(i).border_up ;
   mesh(i).faces(4).zmin = mesh(i).border_south ;
   mesh(i).faces(4).zmax = mesh(i).border_north ;
   mesh(i).faces(4).square_x = [mesh(i).faces(4).xmin ;mesh(i).faces(4).xmax; mesh(i).faces(4).xmax; mesh(i).faces(4).xmin; mesh(i).faces(4).xmin ;];
   mesh(i).faces(4).square_z = [mesh(i).faces(4).zmin ;mesh(i).faces(4).zmin; mesh(i).faces(4).zmax; mesh(i).faces(4).zmax; mesh(i).faces(4).zmin; ];
   mesh(i).faces(4).square_y = [mesh(i).faces(4).ymin ;mesh(i).faces(4).ymax; mesh(i).faces(4).ymax; mesh(i).faces(4).ymin; mesh(i).faces(4).ymin ;];
   
   
    %% 5
   mesh(i).faces(5).xmin = mesh(i).border_left;
   mesh(i).faces(5).xmax = mesh(i).border_left;
   mesh(i).faces(5).ymin = mesh(i).border_down;
   mesh(i).faces(5).ymax = mesh(i).border_up ;
   mesh(i).faces(5).zmin = mesh(i).border_south ;
   mesh(i).faces(5).zmax = mesh(i).border_north ;
   mesh(i).faces(5).square_y = [mesh(i).faces(5).ymin ;mesh(i).faces(5).ymax; mesh(i).faces(5).ymax; mesh(i).faces(5).ymin; mesh(i).faces(5).ymin ;];
   mesh(i).faces(5).square_z = [mesh(i).faces(5).zmin ;mesh(i).faces(5).zmin; mesh(i).faces(5).zmax; mesh(i).faces(5).zmax; mesh(i).faces(5).zmin; ];
   mesh(i).faces(5).square_x = [mesh(i).faces(5).xmin ;mesh(i).faces(5).xmax; mesh(i).faces(5).xmax; mesh(i).faces(5).xmin; mesh(i).faces(5).xmin ;];
   
   
   
   
    %% 6
   mesh(i).faces(6).xmin = mesh(i).border_left;
   mesh(i).faces(6).xmax = mesh(i).border_right;
   mesh(i).faces(6).ymin = mesh(i).border_down;
   mesh(i).faces(6).ymax = mesh(i).border_down ;
   mesh(i).faces(6).zmin = mesh(i).border_south ;
   mesh(i).faces(6).zmax = mesh(i).border_north ;
   mesh(i).faces(6).square_x = [mesh(i).faces(6).xmin ;mesh(i).faces(6).xmax; mesh(i).faces(6).xmax; mesh(i).faces(6).xmin; mesh(i).faces(6).xmin ;];
   mesh(i).faces(6).square_z = [mesh(i).faces(6).zmin ;mesh(i).faces(6).zmin; mesh(i).faces(6).zmax; mesh(i).faces(6).zmax; mesh(i).faces(6).zmin; ];
   mesh(i).faces(6).square_y = [mesh(i).faces(6).ymin ;mesh(i).faces(6).ymax; mesh(i).faces(6).ymax; mesh(i).faces(6).ymin; mesh(i).faces(6).ymin ;];
  
   
end

% close(h);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh = postprocessing_neighbors(mesh)


n_mesh = length(mesh);
for i =1 : n_mesh
    mesh(i).numero = i;
end;


% h       = waitbar(0, 'Post processing neighbors ....');

%% stupid searches for neighbors
for i   = 1 : n_mesh
%     waitbar(i/n_mesh,h);
    mesh(i).neighbors = [];
    
    %% 1
    xxi = mesh(i).faces(1).square_x;
    yyi = mesh(i).faces(1).square_y;
    loc = [];
    
    for j = 1 : n_mesh
        
        xxj = mesh(j).faces(2).square_x;
        yyj = mesh(j).faces(2).square_y;
        
        II  = inpolygon( xxi, yyi, xxj, yyj);
        JJ  = inpolygon( xxj, yyj, xxi, yyi);
        
        if ((sum(II)==5)|(sum(JJ)==5))&(mesh(i).faces(1).square_z==mesh(j).faces(2).square_z)
           loc = [loc,j ];
        end
        
    end
    mesh(i).neighbors = [mesh(i).neighbors ,loc];
    for kkkk = 1 : length(loc)
        mesh(i).plus_z(kkkk) = mesh(loc(kkkk)).center_z;
    end
    mesh(i).plus_z_index = loc;
    
    
    
    %% 2
    xxi = mesh(i).faces(2).square_x;
    yyi = mesh(i).faces(2).square_y;
    loc = [];
    
    for j = 1 : n_mesh
        
        xxj = mesh(j).faces(1).square_x;
        yyj = mesh(j).faces(1).square_y;
        
        II  = inpolygon( xxi, yyi, xxj, yyj);
        JJ  = inpolygon( xxj, yyj, xxi, yyi);
        
        if ((sum(II)==5)|(sum(JJ)==5))&(mesh(i).faces(2).square_z==mesh(j).faces(1).square_z)
           loc = [loc,j ];
        end
        
    end
    mesh(i).neighbors = [mesh(i).neighbors ,loc];
    for kkkk = 1 : length(loc)
        mesh(i).minus_z(kkkk) = mesh(loc(kkkk)).center_z;
    end
    mesh(i).minus_z_index = loc;
    
    
    
    %% 3
    
    yyi = mesh(i).faces(3).square_y;
    zzi = mesh(i).faces(3).square_z;
    loc = [];
    
    for j = 1 : n_mesh
        
        yyj = mesh(j).faces(5).square_y;
        zzj = mesh(j).faces(5).square_z;
        
        II  = inpolygon( yyi, zzi, yyj, zzj);
        JJ  = inpolygon( yyj, zzj, yyi, zzi);
        
        if ((sum(II)==5)|(sum(JJ)==5))&(mesh(i).faces(3).square_x==mesh(j).faces(5).square_x)
           loc = [loc,j ];
        end
        
    end
    mesh(i).neighbors = [mesh(i).neighbors ,loc];
    for kkkk = 1 : length(loc)
        mesh(i).plus_x(kkkk) = mesh(loc(kkkk)).center_x;
    end
    mesh(i).plus_x_index = loc;
    
    
    %% 4
    
    xxi = mesh(i).faces(4).square_x;
    zzi = mesh(i).faces(4).square_z;
    loc = [];
    
    for j = 1 : n_mesh
        
        xxj = mesh(j).faces(6).square_x;
        zzj = mesh(j).faces(6).square_z;
        
        II  = inpolygon( xxi, zzi, xxj, zzj);
        JJ  = inpolygon( xxj, zzj, xxi, zzi);
        
        if ((sum(II)==5)|(sum(JJ)==5))&(mesh(i).faces(4).square_y == mesh(j).faces(6).square_y)
           loc = [loc,j ];
        end
        
    end
    mesh(i).neighbors = [mesh(i).neighbors ,loc];
    for kkkk = 1 : length(loc)
        mesh(i).plus_y(kkkk) = mesh(loc(kkkk)).center_y;
    end
    mesh(i).plus_y_index = loc;
    
    
    %% 5
    
    yyi = mesh(i).faces(5).square_y;
    zzi = mesh(i).faces(5).square_z;
    loc = [];
    
    for j = 1 : n_mesh
        
        yyj = mesh(j).faces(3).square_y;
        zzj = mesh(j).faces(3).square_z;
        
        II  = inpolygon( yyi, zzi, yyj, zzj);
        JJ  = inpolygon( yyj, zzj, yyi, zzi);
        
        if ((sum(II)==5)|(sum(JJ)==5))&(mesh(i).faces(5).square_x==mesh(j).faces(3).square_x)
           loc = [loc,j ];
        end
        
    end
    mesh(i).neighbors = [mesh(i).neighbors ,loc];
    for kkkk = 1 : length(loc)
        mesh(i).minus_x(kkkk) = mesh(loc(kkkk)).center_x;
    end
    mesh(i).minus_x_index = loc;
    
    %% 6
    
    xxi = mesh(i).faces(6).square_x;
    zzi = mesh(i).faces(6).square_z;
    loc = [];
    
    for j = 1 : n_mesh
        
        xxj = mesh(j).faces(4).square_x;
        zzj = mesh(j).faces(4).square_z;
        
        II  = inpolygon( xxi, zzi, xxj, zzj);
        JJ  = inpolygon( xxj, zzj, xxi, zzi);
        
        if ((sum(II)==5)|(sum(JJ)==5))&(mesh(i).faces(6).square_y==mesh(j).faces(4).square_y)
           loc = [loc,j ];
        end
        
    end
    mesh(i).neighbors = [mesh(i).neighbors ,loc];
    for kkkk = 1 : length(loc)
        mesh(i).minus_y(kkkk) = mesh(loc(kkkk)).center_y;
    end
    mesh(i).minus_y_index = loc;
    %%
    
    

end

% close(h);

for i =1 : n_mesh
   
    mesh(i).neighbors = unique(mesh(i).neighbors);
    
end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_tot, y_tot,z_tot, t_tot,dx_tot , dy_tot, dz_tot] = build_point_map(tout)

x_tot  = tout(:,1);
y_tot  = tout(:,2);
z_tot  = tout(:,3);

t_tot  = tout(:,4);

dx_tot = tout(:,5);
dy_tot = tout(:,6);
dz_tot = tout(:,7);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mesh, movie_per_frame] = postprocessing_trajectories(mesh, x, y, t, movie_per_frame )


for i = 1 : length(mesh)
     
    II      = find(       ( x >= mesh(i).border_left  )...
                        & ( x <= mesh(i).border_right )...
                        & ( y <= mesh(i).border_up    )...
                        & ( y >= mesh(i).border_down  )...
                        & ( z <= mesh(i).border_north    )...
                        & ( z >= mesh(i).border_south  ) );

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