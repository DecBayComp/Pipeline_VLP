function [Maps] = give_D_f(parameters,  Maps,dt)



if parameters.d == 2

    
    n_zone                          = length(Maps);
    sigma                           = parameters.sigma;
    opt                             = parameters.opt;
    param_out                       = zeros(n_zone,3);
    Feval                           = zeros(n_zone,1);

    for i = 1 : n_zone  
        % fprintf('optimizing D, fx, fy in %i\n', i);
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



elseif parameters.d == 3
    
    n_zone                          = length(Maps);
    sigma                           = parameters.sigma;
    opt                             = parameters.opt;
    param_out                       = zeros(n_zone,4);
    Feval                           = zeros(n_zone,1);

    for i = 1 : n_zone  
        % fprintf('optimizing D, fx, fy, fz in %i\n', i);
        dx    = Maps(i).dx;
        dy    = Maps(i).dy;
        dz    = Maps(i).dz;

        dr2   = dx.^2 + dy.^2 + dz.^2;

%         if isempty(dr2)
        if length(dx)<3
            param_out(i,:) = NaN(1,4);
            Feval(i) = -inf;
        else

        nn                  = length(dx);
        D_init              = sum(dr2./(6*dt))./nn;
        fx_init             = sum(dx./(D_init*dt))./nn;
        fy_init             = sum(dy./(D_init*dt))./nn;
        fz_init             = sum(dz./(D_init*dt))./nn;

        param_init          = [D_init, fx_init, fy_init, fz_init];
        [param_out(i,:), Feval(i)]  = fminunc(@(param) likelihood_D_f_3D(dx,dy,dz, dt, param,sigma),param_init, opt); 
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




end

        
  

end