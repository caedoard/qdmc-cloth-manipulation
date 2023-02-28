function [phi,dphi] = simulate_cloth_step(x,u,model)

szmdl = length(x)/2;
phi1  = x(1:szmdl);
dphi1 = x(szmdl+1:2*szmdl);

dt = model.dt;

[phi,dphi] = implicit_euler(dt, phi1, dphi1, u, ...
                            model.Fg, model.rho, model.A_b, ...
                            model.Minv, model.D, model.C, ...
                            model.dE, model.Mlum, ...
                            model.n_conds, model.n_nodos, ...
                            model.Cphi0, model.den, model.bnd_aux);

end