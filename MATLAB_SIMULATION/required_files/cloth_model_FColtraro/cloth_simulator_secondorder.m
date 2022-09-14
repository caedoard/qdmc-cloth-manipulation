function [xdot] = cloth_simulator_secondorder(x,u,model_4x4)


ll = length(x);
phi1 = x(1:ll/4);
dphi1 = x((1+ll/4):ll/2);

phi0 = x((1+ll/2):0.75*ll);
dphi0 = x((1+0.75*ll):ll);

            
[phi,dphi] = imex_franco(model_4x4.dt,phi0,dphi0,phi1,dphi1,...
                           model_4x4.Fg,model_4x4.rho,model_4x4.A_b,model_4x4.Minv,model_4x4.Mcons,model_4x4.K,model_4x4.D,model_4x4.C,...
                           model_4x4.n_conds,model_4x4.n_nodos,model_4x4.Cphi0,u);



xdot=[phi;dphi;phi1;dphi1];

end

