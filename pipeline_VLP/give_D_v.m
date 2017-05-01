function [Maps] = give_D_v(parameters,  Maps,dt)

% selective_kill_parallel_multi_version;
% start_check_parallel_multi_version;


        
if parameters.d == 2
    
    
    n_zone                          = length(Maps);
    sigma                           = parameters.sigma;
    opt                             = parameters.opt;
    param_out                       = zeros(n_zone,3);
    Feval                           = zeros(n_zone,1);

    for i = 1 : n_zone  
       %  fprintf('optimizing D, vx, vy in %i\t length x %i\n', i, length(Maps(i).dx));
        dx    = Maps(i).dx;
        dy    = Maps(i).dy;
        dr2   = dx.^2 + dy.^2;

%         if isempty(dr2)
        if length(dx)<3
            param_out(i,:) = NaN(1,3);
            Feval(i) = -inf;
        else

            nn    = length(dx);
            D_init              = sum(dr2./(4*dt))./nn;
            vx_init             = sum(dx./(dt))./nn;
            vy_init             = sum(dy./(dt))./nn;
            param_init          = [D_init, vx_init, vy_init];
            [param_out(i,:), Feval(i)]  = fminunc(@(param) likelihood_D_fb(dx,dy, dt, param,sigma),param_init, opt); 
        end
    end
%     fprintf('here 1\n');
    for i = 1 : n_zone

        Maps(i).D     = param_out(i,1);
        Maps(i).vx    = param_out(i,2);
        Maps(i).vy    = param_out(i,3);
        Maps(i).v     = sqrt(Maps(i).vx.^2 + Maps(i).vy.^2);
        Maps(i).log_likelihood     = -Feval(i);


    end
%     fprintf('here 2\n');
elseif parameters.d ==3
    
     n_zone                          = length(Maps);
    sigma                           = parameters.sigma;
    opt                             = parameters.opt;
    param_out                       = zeros(n_zone,4);
    Feval                           = zeros(n_zone,1);

    for i = 1 : n_zone  
        % fprintf('optimizing D, vx, vy, vz in %i\n', i);
        dx    = Maps(i).dx;
        dy    = Maps(i).dy;
        dz    = Maps(i).dz;

        dr2   = dx.^2 + dy.^2 + dz.^2;

%         if isempty(dr2)
        if length(dx)<3
            param_out(i,:) = NaN(1,4);
            Feval(i) = -inf;
        else

            nn    = length(dx);
            D_init              = sum(dr2./(6*dt))./nn;
            vx_init             = sum(dx./(dt))./nn;
            vy_init             = sum(dy./(dt))./nn;
            vz_init             = sum(dz./(dt))./nn;

            param_init          = [D_init, vx_init, vy_init, vz_init];
            [param_out(i,:), Feval(i)]  = fminunc(@(param) likelihood_D_fb_3D(dx,dy,dz, dt, param,sigma),param_init, opt); 
        end
    end

    for i = 1 : n_zone

        Maps(i).D     = param_out(i,1);
        Maps(i).vx    = param_out(i,2);
        Maps(i).vy    = param_out(i,3);
        Maps(i).vz    = param_out(i,4);

        Maps(i).v     = sqrt(Maps(i).vx.^2 + Maps(i).vy.^2 +  Maps(i).vz.^2 );
        Maps(i).log_likelihood     = -Feval(i);


    end


end
        
        
  

end