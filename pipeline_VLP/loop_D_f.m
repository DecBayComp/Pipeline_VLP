function [Maps] = loop_D_f(parameters, tout, Maps,dt)

% kill_parallel_multi_version;
% start_check_parallel_multi_version;

indice = 0; 

  %  while (1)
  if parameters.d == 2
      
      while(indice < parameters.number_stop_loop_tree )
        n_zone                          = length(Maps);

        sigma                           = parameters.sigma;
        opt                             = parameters.opt;
        param_out                       = zeros(n_zone,3);
        Feval                           = zeros(n_zone,1);

        for i = 1 : n_zone  
            % fprintf('optimizing D, fx, fy  in %i\n', i);
            dx    = Maps(i).dx;
            dy    = Maps(i).dy;
            dr2   = dx.^2 + dy.^2;
            nn    = length(dx);
            
            
%         if isempty(dr2)
            if length(dx)<3
                
                param_out(i,:) = [parameters.D_init, NaN, NaN];
                Feval(i) = -inf;
                
            else
                
                D_init              = sum(dr2./(4*dt))./nn;
                fx_init             = sum(dx./(D_init*dt))./nn;
                fy_init             = sum(dy./(D_init*dt))./nn;
                param_init          = [D_init, fx_init, fy_init];
                [param_out(i,:), Feval(i)]  = fminunc(@(param) likelihood_D_f(dx,dy, dt, param,sigma),param_init, opt); 
                
                
            end
            
        end
    
        for i = 1 : n_zone
                
            Maps(i).D     = param_out(i,1);
            Maps(i).fx    = param_out(i,2);
            Maps(i).fy    = param_out(i,3);
            Maps(i).f     = sqrt(Maps(i).fx.^2 + Maps(i).fy.^2);
            Maps(i).log_likelihood     = -Feval(i);
                
                
        end
        
        
        
        Maps        = switch_status_state(Maps,1);
        [Map_points_new]  = build_tree(tout, parameters, Maps);
        Maps        = switch_status_state(Maps, 0);

    
        if (isequaln(Map_points_new, Maps))
            break;
        else
            Maps = Map_points_new;
        end
        
        indice = indice +1 ;
      end
    
  elseif parameters.d == 3
      
      while(indice < parameters.number_stop_loop_tree )
        n_zone                          = length(Maps);

        sigma                           = parameters.sigma;
        opt                             = parameters.opt;
        param_out                       = zeros(n_zone,4);
        Feval                           = zeros(n_zone,1);

        for i = 1 : n_zone  
            % fprintf('optimizing D, fx, fy  in %i\n', i);
            dx    = Maps(i).dx;
            dy    = Maps(i).dy;
            dz    = Maps(i).dz;
            
            dr2   = dx.^2 + dy.^2 + dz.^2;
            nn    = length(dx);
            
            
%         if isempty(dr2)
            if length(dx)<3
                
                param_out(i,:) = [parameters.D_init, NaN, NaN, NaN];
                Feval(i) = -inf;
                
            else
                
                D_init              = sum(dr2./(4*dt))./nn;
                fx_init             = sum(dx./(D_init*dt))./nn;
                fy_init             = sum(dy./(D_init*dt))./nn;
                fz_init             = sum(dz./(D_init*dt))./nn;
                
                param_init          = [D_init, fx_init, fy_init, fz_init];
                [param_out(i,:), Feval(i)]  = fminunc(@(param) likelihood_D_f_3D(dx,dy, dz, dt, param,sigma),param_init, opt); 
                
                
            end
            
        end
    
        for i = 1 : n_zone
                
            Maps(i).D     = param_out(i,1);
            Maps(i).fx    = param_out(i,2);
            Maps(i).fy    = param_out(i,3);
            Maps(i).fz    = param_out(i,4);
            
            Maps(i).f     = sqrt(Maps(i).fx.^2 + Maps(i).fy.^2 +  Maps(i).fz.^2);
            Maps(i).log_likelihood     = -Feval(i);
                
                
        end
        
        
        
        Maps        = switch_status_state(Maps,1);
        [Map_points_new]  = build_tree_3D(tout, parameters, Maps);
        Maps        = switch_status_state(Maps, 0);

    
        if (isequaln(Map_points_new, Maps))
            break;
        else
            Maps = Map_points_new;
        end
        
        indice = indice +1 ;
      end
      
      
  end
    










end