function [ Map_points2 ] = smooth_maps_multi_methods( Maps, field , method, scale, varargin )

%% check if multiple processors can be used

% selective_kill_parallel_multi_version;
% start_check_parallel_multi_version;
    
if (nargin <= 3)
    scale = 0.2;
end
    
 Map_points2 =  Maps ;
 
if ~isfield(Maps,'neighbors');
else 
 
 switch lower(method)
     
     case 'neighbors'
          
         for i =1 : length(Maps)
 
 
            switch lower(field)
         
                case 'diffusion'
             
                D_loc = [Maps(i).D];             
                for j = 1 : length( Maps(i).neighbors )
                    indice = Maps(i).neighbors(j);
                    if (isnan(Maps(indice).D))
                    else
                        D_loc =  [D_loc ;Maps(indice).D ];
                    end
                 
                end
                Map_points2(i).D = mean(D_loc);
             
                case 'potential'
             
                V_loc = [Maps(i).V];             
                for j = 1 : length( Maps(i).neighbors )
                    indice = Maps(i).neighbors(j);
                    if (isnan(Maps(indice).D))
                    else
                        V_loc =  [V_loc ;Maps(indice).V ];
                    end
             
                end
                Map_points2(i).V = mean(V_loc);
             
                case 'all'
             
                D_loc = [Maps(i).D];
                V_loc = [Maps(i).V];
                f_loc = [Maps(i).f];
                fx_loc =  [Maps(i).fx];
                fy_loc =  [Maps(i).fy];  
                v_loc = [Maps(i).v];
                vx_loc =  [Maps(i).vx];
                vy_loc =  [Maps(i).vy]; 
                
                for j = 1 : length( Maps(i).neighbors )
                    indice = Maps(i).neighbors(j);
                
                    if (isnan(Maps(indice).D))
                    else
                        D_loc =  [D_loc ;Maps(indice).D ];
                    end
                    
                    if (isnan(Maps(indice).D))
                    else
                        V_loc =  [V_loc ;Maps(indice).V ];
                    end
                    
%                     V_loc =  [V_loc ;Map_points(indice).V ];
%                     D_loc =  [D_loc ;Map_points(indice).D ];
                    if (isnan(Maps(indice).f))
                    else
                        f_loc =  [f_loc ;Maps(indice).f ];
                        fx_loc =  [fx_loc ;Maps(indice).fx ];
                        fy_loc =  [fy_loc ;Maps(indice).fy ];
                    end
                    if (isnan(Maps(indice).v))
                    else
                        v_loc =  [v_loc ;Maps(indice).v ];
                        vx_loc =  [vx_loc ;Maps(indice).vx ];
                        vy_loc =  [vy_loc ;Maps(indice).vy ];
                    end
                 
                end
                Map_points2(i).D = mean(D_loc);
                Map_points2(i).V = mean(V_loc);
                Map_points2(i).f = mean(f_loc); 
                Map_points2(i).fx = mean(fx_loc); 
                Map_points2(i).fy = mean(fy_loc);
                Map_points2(i).v = mean(v_loc); 
                Map_points2(i).vx = mean(vx_loc); 
                Map_points2(i).vy = mean(vy_loc); 
  
                case 'forces'
             
                f_loc = [Maps(i).f];
                fx_loc =  [Maps(i).fx];
                fy_loc =  [Maps(i).fy];          
                for j = 1 : length( Maps(i).neighbors )
                    indice = Maps(i).neighbors(j);
                
                    if (isnan(Maps(indice).f))
                    else
                        f_loc =  [f_loc ;Maps(indice).f ];
                        fx_loc =  [fx_loc ;Maps(indice).fx ];
                        fy_loc =  [fy_loc ;Maps(indice).fy ];
                    end
                
                end
                Map_points2(i).f = mean(f_loc); 
                Map_points2(i).fx = mean(fx_loc); 
                Map_points2(i).fy = mean(fy_loc); 
                
                case 'drifts'
             
                v_loc = [Maps(i).v];
                vx_loc =  [Maps(i).vx];
                vy_loc =  [Maps(i).vy];          
                for j = 1 : length( Maps(i).neighbors )
                    indice = Maps(i).neighbors(j);
                
                    if (isnan(Maps(indice).v))
                    else
                        v_loc =  [v_loc ;Maps(indice).v ];
                        vx_loc =  [vx_loc ;Maps(indice).vx ];
                        vy_loc =  [vy_loc ;Maps(indice).vy ];
                    end
                
                end
                Map_points2(i).v = mean(v_loc); 
                Map_points2(i).vx = mean(vx_loc); 
                Map_points2(i).vy = mean(vy_loc); 
                
                
               
            end
     
        end
         
         
         
     case 'gaussian'
         
         
         
         
     %    for i =1 : length(Map_points)
            
 
            switch lower(field)
         
                case 'diffusion'
                    
                    
                    for i = 1 : length(Maps)
                        
                        x(i,1) = Maps(i).center_x;
                        y(i,1) = Maps(i).center_y;
                        D(i,1) = Maps(i).D;

                    end
                    
                    for i = 1 : length(Maps)
                        
                        r2                = ( (x - x(i,1)).^2 + (y - y(i,1)).^2 );
                        D_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*D )./nansum(exp(-r2./(2.*scale.^2)).*D./D);
                        Map_points2(i).D = D_eff(i,1);
                        
                        
                    end
                
             
                case 'potential'
             
                    for i = 1 : length(Maps)
                        
                        x(i,1) = Maps(i).center_x;
                        y(i,1) = Maps(i).center_y;
                        V(i,1) = Maps(i).V;

                    end
                    
                    for i = 1 : length(Maps)
                        
                        r2                = ( (x - x(i,1)).^2 + (y - y(i,1)).^2 );
                        V_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*V )./nansum(exp(-r2./(2.*scale.^2)).*V./V);
                        Map_points2(i).V = V_eff(i,1);
                        
                    end
                    
                
             
                case 'all'
                    
                    for i = 1 : length(Maps)
                        
                        x(i,1) = Maps(i).center_x;
                        y(i,1) = Maps(i).center_y;
                        f(i,1) = Maps(i).f;
                        fx(i,1) = Maps(i).fx;
                        fy(i,1) = Maps(i).fy;
                        v(i,1) = Maps(i).v;
                        vx(i,1) = Maps(i).vx;
                        vy(i,1) = Maps(i).vy;
                        V(i,1) = Maps(i).V;
                        D(i,1) = Maps(i).D;
                        
                    end
                    
                    for i = 1 : length(Maps)
                        
                        r2                = ( (x - x(i,1)).^2 + (y - y(i,1)).^2 );
                        f_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*f )./nansum(exp(-r2./(2.*scale.^2)).*f./f);
                        fx_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*fx )./nansum(exp(-r2./(2.*scale.^2)).*fx./fx);
                        fy_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*fy )./nansum(exp(-r2./(2.*scale.^2)).*fy./fy);
                        v_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*v )./nansum(exp(-r2./(2.*scale.^2)).*v./v);
                        vx_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*vx )./nansum(exp(-r2./(2.*scale.^2)).*vx./vx);
                        vy_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*vy )./nansum(exp(-r2./(2.*scale.^2)).*vy./vy);
                        V_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*V )./nansum(exp(-r2./(2.*scale.^2)).*V./V);
                        D_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*D )./nansum(exp(-r2./(2.*scale.^2)).*D./D);
                         
                        Map_points2(i).f = f_eff(i,1);
                        Map_points2(i).fx = fx_eff(i,1);
                        Map_points2(i).fy = fy_eff(i,1);
                        Map_points2(i).v = v_eff(i,1);
                        Map_points2(i).vx = vx_eff(i,1);
                        Map_points2(i).vy = vy_eff(i,1);
                        Map_points2(i).V = V_eff(i,1);
                        Map_points2(i).D = D_eff(i,1);
                        
                    end
                    
             
                case 'forces'
             
                    for i = 1 : length(Maps)
                        
                        x(i,1) = Maps(i).center_x;
                        y(i,1) = Maps(i).center_y;
                        f(i,1) = Maps(i).f;
                        fx(i,1) = Maps(i).fx;
                        fy(i,1) = Maps(i).fy;
                        
                    end
                    
                    for i = 1 : length(Maps)
                        
                        r2                = ( (x - x(i,1)).^2 + (y - y(i,1)).^2 );
                        
                        f_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*f )./nansum(exp(-r2./(2.*scale.^2)).*f./f);
                        fx_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*fx )./nansum(exp(-r2./(2.*scale.^2)).*fx./fx);
                        fy_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*fy )./nansum(exp(-r2./(2.*scale.^2)).*fy./fy);
                        
                        Map_points2(i).f = f_eff(i,1);
                        Map_points2(i).fx = fx_eff(i,1);
                        Map_points2(i).fy = fy_eff(i,1);
                        
                        
                        
                    end

                    
                case 'drifts'
             
                     for i = 1 : length(Maps)
                        
                        x(i,1) = Maps(i).center_x;
                        y(i,1) = Maps(i).center_y;
                        v(i,1) = Maps(i).v;
                        vx(i,1) = Maps(i).vx;
                        vy(i,1) = Maps(i).vy;
                        
                    end
                    
                    for i = 1 : length(Maps)
                        
                        r2                = ( (x - x(i,1)).^2 + (y - y(i,1)).^2 );
                        
                        v_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*v )./nansum(exp(-r2./(2.*scale.^2)).*v./v);
                        vx_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*vx )./nansum(exp(-r2./(2.*scale.^2)).*vx./vx);
                        vy_eff(i,1)       = nansum(exp(-r2./(2.*scale.^2)).*vy )./nansum(exp(-r2./(2.*scale.^2)).*vy./vy);
                        
                        Map_points2(i).v = v_eff(i,1);
                        Map_points2(i).vx = vx_eff(i,1);
                        Map_points2(i).vy = vy_eff(i,1);
                        
                        
                        
                    end
                    
                     
               
            end
     
    %    end
         
         
         
         
         
     case 'radius'
         
         
         switch lower(field)
         
                case 'diffusion'
                    
                    
                    for i = 1 : length(Maps)
                        
                        x(i,1) = Maps(i).center_x;
                        y(i,1) = Maps(i).center_y;
                        D(i,1) = Maps(i).D;

                    end
                    
                    for i = 1 : length(Maps)
                        
                        r2                = ( (x - x(i,1)).^2 + (y - y(i,1)).^2 );
                        II                = r2 <= scale^2;
                        D_eff(i,1)       = nanmean(D(II));
                        Map_points2(i).D = D_eff(i,1);
                        
                        
                    end
                
             
                case 'potential'
             
                    for i = 1 : length(Maps)
                        
                        x(i,1) = Maps(i).center_x;
                        y(i,1) = Maps(i).center_y;
                        V(i,1) = Maps(i).V;

                    end
                    
                    for i = 1 : length(Maps)
                        
                        r2                = ( (x - x(i,1)).^2 + (y - y(i,1)).^2 );
                        II                = r2 <= scale^2;
                        V_eff(i,1)       = nanmean(V(II));
                        Map_points2(i).V = V_eff(i,1);
                        
                    end
                    
                
             
                case 'all'
                    
                    for i = 1 : length(Maps)
                        
                        x(i,1) = Maps(i).center_x;
                        y(i,1) = Maps(i).center_y;
                        f(i,1) = Maps(i).f;
                        fx(i,1) = Maps(i).fx;
                        fy(i,1) = Maps(i).fy;
                        v(i,1) = Maps(i).v;
                        vx(i,1) = Maps(i).vx;
                        vy(i,1) = Maps(i).vy;
                        V(i,1) = Maps(i).V;
                        D(i,1) = Maps(i).D;
                        
                    end
                    
                    for i = 1 : length(Maps)
                        
                        r2                = ( (x - x(i,1)).^2 + (y - y(i,1)).^2 );
                        II                = r2 <= scale^2;
                        f_eff(i,1)        = nanmean(f(II));
                        fx_eff(i,1)       = nanmean(fx(II));
                        fy_eff(i,1)       = nanmean(fy(II));
                        v_eff(i,1)        = nanmean(v(II));
                        vx_eff(i,1)       = nanmean(vx(II));
                        vy_eff(i,1)       = nanmean(vy(II));
                        V_eff(i,1)        = nanmean(V(II));
                        D_eff(i,1)        = nanmean(D(II));
                         
                        Map_points2(i).f = f_eff(i,1);
                        Map_points2(i).fx = fx_eff(i,1);
                        Map_points2(i).fy = fy_eff(i,1);
                        Map_points2(i).v = v_eff(i,1);
                        Map_points2(i).vx = vx_eff(i,1);
                        Map_points2(i).vy = vy_eff(i,1);
                        Map_points2(i).V = V_eff(i,1);
                        Map_points2(i).D = D_eff(i,1);
                        
                    end
                    
             
                
  
                case 'forces'
             
                     for i = 1 : length(Maps)
                        
                        x(i,1) = Maps(i).center_x;
                        y(i,1) = Maps(i).center_y;
                        f(i,1) = Maps(i).f;
                        fx(i,1) = Maps(i).fx;
                        fy(i,1) = Maps(i).fy;
                        
                    end
                    
                    for i = 1 : length(Maps)
                        
                        r2                = ( (x - x(i,1)).^2 + (y - y(i,1)).^2 );
                        II                = r2 <= scale^2;
                        f_eff(i,1)        = nanmean(f(II));
                        fx_eff(i,1)       = nanmean(fx(II));
                        fy_eff(i,1)       = nanmean(fy(II));

                        Map_points2(i).f = f_eff(i,1);
                        Map_points2(i).fx = fx_eff(i,1);
                        Map_points2(i).fy = fy_eff(i,1);
                        
                        
                        
                    end
                    
                
                    
                case 'drifts'
             
                     for i = 1 : length(Maps)
                        
                        x(i,1) = Maps(i).center_x;
                        y(i,1) = Maps(i).center_y;
                        v(i,1) = Maps(i).v;
                        vx(i,1) = Maps(i).vx;
                        vy(i,1) = Maps(i).vy;
                        
                    end
                    
                    for i = 1 : length(Maps)
                        
                        r2                = ( (x - x(i,1)).^2 + (y - y(i,1)).^2 );
                        II                = r2 <= scale^2;
                        v_eff(i,1)        = nanmean(v(II));
                        vx_eff(i,1)       = nanmean(vx(II));
                        vy_eff(i,1)       = nanmean(vy(II));

                        Map_points2(i).v = v_eff(i,1);
                        Map_points2(i).vx = vx_eff(i,1);
                        Map_points2(i).vy = vy_eff(i,1);
                        
                        
                        
                    end                         
                    
                    
                    
               
            end
         
         
         
         
 end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
