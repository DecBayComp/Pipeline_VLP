function [Maps] = loop_D(parameters, tout, Maps,dt)



indice = 0; 

    if parameters.d == 2
        
        while(indice < parameters.number_stop_loop_tree )
            n_zone                          = length(Maps);
            sigma                           = parameters.sigma;
            opt                             = parameters.opt;
            D_out                           = zeros(n_zone,1);
            Feval                           = zeros(n_zone,1);
        
            for i = 1 : n_zone  
                % fprintf('optimizing D in %i\n', i);
                dr2    = Maps(i).dx.^2 + Maps(i).dy.^2;
                D_init = sum(dr2/(4.*dt))./length(dr2);
                Dmax    = 2*D_init;
                Dmin    = D_init./5;
            
                if isempty(dr2);
                    D_out(i) = D_init;
                    Feval(i) = -inf;
                else
                    [D_out(i), Feval(i)]                 = fminbnd(@(D) likelihood_D(dr2,dt, D,sigma),Dmin, Dmax, opt);            
                end

            end
    
            for i =1 : n_zone
                Maps(i).D                = D_out(i);
                Maps(i).log_likelihood      = -Feval(i);
            end
        
            Maps        = switch_status_state(Maps,1);
            [Map_points_new]  = build_tree(tout, parameters, Maps);
            Maps        = switch_status_state(Maps, 0);
    
            if (isequaln(Map_points_new, Maps))
                break;
            else
                Maps = Map_points_new;
            end
            %isequaln(quadtree_mesh, quadtree_mesh);
            indice = indice +1 ;
        end
        
    elseif parameters.d ==3
        
        while(indice < parameters.number_stop_loop_tree )
            n_zone                          = length(Maps);
            sigma                           = parameters.sigma;
            opt                             = parameters.opt;
            D_out                           = zeros(n_zone,1);
            Feval                           = zeros(n_zone,1);
        
            for i = 1 : n_zone  
                % fprintf('optimizing D in %i\n', i);
                dr2    = Maps(i).dx.^2 + Maps(i).dy.^2 + + Maps(i).dz.^2;
                D_init = sum(dr2/(6.*dt))./length(dr2);
                Dmax    = 2*D_init;
                Dmin    = D_init./5;
            
                if isempty(dr2);
                    D_out(i) = D_init;
                    Feval(i) = -inf;
                else
                    [D_out(i), Feval(i)]                 = fminbnd(@(D) likelihood_D_3D(dr2,dt, D,sigma),Dmin, Dmax, opt);            
                end

            end
    
            for i =1 : n_zone
                Maps(i).D                = D_out(i);
                Maps(i).log_likelihood      = -Feval(i);
            end
        
            Maps        = switch_status_state(Maps,1);
            [Map_points_new]  = build_tree_3D(tout, parameters, Maps);
            Maps        = switch_status_state(Maps, 0);
    
            if (isequaln(Map_points_new, Maps))
                break;
            else
                Maps = Map_points_new;
            end
            %isequaln(quadtree_mesh, quadtree_mesh);
            indice = indice +1 ;
        end
        
        
    end

    


end