function [Cnew,n_conds] = processConstraints(C,n_nodos,nodes_int,esquinas)
Cu = C.Cu; Cv = C.Cv; Cuv = C.Cuv;
%solo nos quedamos con las condiciones interiores y la condicion para
%preservar el angulo de la esquina
Cnew = [Cu(:,nodes_int),Cv(:,nodes_int),Cuv(:,[nodes_int,esquinas])]; 
n_conds = size(Cnew,2);
Cnew = reshape(Cnew',[n_conds*n_nodos,n_nodos]);
end


