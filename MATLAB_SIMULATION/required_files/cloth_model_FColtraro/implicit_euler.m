function [phi,dphi,cost,cost1,n_iter] = implicit_euler(dt,phi1,dphi1,u_t,...
                                                       Fg,rho,A_b,Minv,D,C,dE,Mlum,...
                                                       n_conds,n_nodos,C0,den,bnd_aux)
%paso sin restricciones
cte = (dt^2/rho); 
Fimp = (dt^2)*Fg + dt*rho*Mlum*dphi1 + (rho*Mlum + dt*D)*phi1;
phi = dE\Fimp;
%restriccion
[Cphi,J,tmp] = fun_C(phi,C,A_b,bnd_aux,n_nodos,n_conds);
C0((end-(length(u_t)-1)):end) = u_t; error = Cphi-C0; 
cost = tmp; cost1 = 0;

%iteraciones
error_max = 2;
n_iter = 0; 
t0_iter = tic;
t_iter = toc(t0_iter);
timeout = 0.1;
while (error_max > 0.1 || n_iter < 1) && t_iter < timeout
    MinvJt = cte*Minv*J'; A = J*MinvJt;
    A = 0.5*(A+A');
    t1 = tic;
    dlt_lambda = A\error;
    cost1 = cost1 + toc(t1);
    dlt_phi = -MinvJt*dlt_lambda;
    phi = phi + dlt_phi;
    [Cphi,J,tmp] = fun_C(phi,C,A_b,bnd_aux,n_nodos,n_conds);
    error = Cphi-C0;
    cost = cost + tmp;
    error_max = 100*max(abs(error./den));
    n_iter = n_iter+1;
    t_iter = toc(t0_iter);
end
%actualizamos la velocidad
dphi = (phi-phi1)/dt;

end
