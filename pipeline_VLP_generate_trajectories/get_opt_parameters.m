function parameters = get_opt_parameters(parameters)

%%
try
parameters.opt = optimset('Display','none', ... % 'off','final','notify'
'MaxFunEvals',1e7,'MaxIter',1e7, ...
'TolX',1e-7,'TolFun',1e-7,'LargeScale','off');
end
%%
try
parameters.opt_loc = optimset('Display','iter', ... % 'off','final','notify'
'MaxFunEvals',1e7,'MaxIter',1e7, ...
'TolX',1e-7,'TolFun',1e-7,'LargeScale','off');
end
%%
try
parameters.opt2 = optimset('Display','iter', ... % 'off','final','notify'
'MaxFunEvals',1e7,'MaxIter',1e7, ...
'TolX',1e-7,'TolFun',1e-7,'LargeScale','off');
end
%%
try
parameters.opts           = statset('MaxIter' ,50, 'Display', 'off');
end
try
parameters.opts_boosting  = statset('MaxIter' ,2, 'Display', 'off');
end
%%
try
parameters.opt3 = psoptimset('TolMesh',1e-6,'TolX',1e-6, 'TolFun', 1e-6,'MaxIter', 1e6, 'MaxFunEvals', 1e9,...
'Display', 'Iter','UseParallel', 'always'); 
end
%%
try
parameters.opt4 = optimset('TolX',1e-6, 'TolFun', 1e-6,'MaxIter', 1e6, 'MaxFunEvals', 1e9,...
'Display', 'Iter'); 
end
%%
 

end