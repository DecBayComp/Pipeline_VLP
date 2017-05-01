%function [x_tot, y_tot, t_tot,dx_tot , dy_tot, z_tot, dz_tot, varargout] = build_point_map(tout)
function [varargout] = build_point_map(tout)
%% extract position, time and displacement from tout
%to be found in :  Matlab/projet/Mapping_without_tracking/prepare_data


[~,m] = size(tout);
if m==5
   
    x_tot       = tout(:,1);
    y_tot       = tout(:,2);
    t_tot       = tout(:,3);
    dx_tot      = tout(:,4);
    dy_tot      = tout(:,5);
    [x_tot, I]  = sort(x_tot);
    y_tot       = y_tot(I);
    t_tot       = t_tot(I);
    dx_tot      = dx_tot(I);
    dy_tot      = dy_tot(I);
    

    varargout(1) = {x_tot};
    varargout(2) = {y_tot};
    varargout(3) = {t_tot};
    varargout(4) = {dx_tot};
    varargout(5) = {dy_tot};
    
else
    
    x_tot       = tout(:,1);
    y_tot       = tout(:,2);
    z_tot       = tout(:,3);
    t_tot       = tout(:,4);
    dx_tot      = tout(:,5);
    dy_tot      = tout(:,6);
    dz_tot      = tout(:,7);

    [x_tot, I]  = sort(x_tot);
    y_tot       = y_tot(I);
    z_tot       = z_tot(I);
    t_tot       = t_tot(I);
    dx_tot      = dx_tot(I);
    dy_tot      = dy_tot(I);
    dz_tot      = dz_tot(I);
    
    varargout(1) = {x_tot};
    varargout(2) = {y_tot};
    varargout(3) = {z_tot};
    varargout(4) = {t_tot};
    varargout(5) = {dx_tot};
    varargout(6) = {dy_tot};    
    varargout(7) = {dz_tot};
    
    

end









end