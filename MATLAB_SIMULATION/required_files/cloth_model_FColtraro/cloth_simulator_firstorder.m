function [xdot] = cloth_simulator_firstorder(x,u,model_4x4)


ll = length(x);
phi0 = x(1:ll/2);
dphi0 = x((1+ll/2):ll);



%para iniciar el integrador de segundo orden
[phi1,dphi1] = explicit_grispun(model_4x4.dt,phi0,dphi0,model_4x4.Fg,model_4x4.rho,...
                        model_4x4.A_b,model_4x4.Mlum,model_4x4.Minv,model_4x4.Mcons,...
                        model_4x4.K,model_4x4.D,model_4x4.C,model_4x4.n_conds,...
                        model_4x4.n_nodos, model_4x4.Cphi0,u);
            


xdot=[phi1;dphi1];

end

