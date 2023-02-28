function [Cphi,DCphi,tmp] = fun_C(phi,C,A_b,bnd_aux,n_nodos,n_conds)
%devuelve la funcion de restricciones y su gradiente 
%evaluada en la superficie actual
Eb = bnd_aux.Eb; n_aristas = size(Eb,1);
t0 = tic;
phi_xyz = reshape(phi,[n_nodos,3]);
%borde
vecs = phi_xyz(Eb(:,2),:) - phi_xyz(Eb(:,1),:); 
longs = dot(vecs,vecs,2);
grad_f = 2*vecs; grad_0 = -2*vecs;
K = [grad_f(:,1)',grad_f(:,2)',grad_f(:,3)',...
     grad_0(:,1)',grad_0(:,2)',grad_0(:,3)'];
grad_bnd = sparse(bnd_aux.I,bnd_aux.J,K,n_aristas,3*n_nodos);
%interior
grad = reshape(C*phi_xyz,[n_conds,3*n_nodos]);
%Valor de la funcion
grad_phi = grad*phi; Aphi = A_b*phi;
Cphi = vertcat(grad_phi,longs,Aphi);
%Jacobiano
DCphi = vertcat(2*grad,grad_bnd,A_b);
tmp = toc(t0);

