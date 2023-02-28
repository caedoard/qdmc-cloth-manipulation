function [C,Mlum,Minv,D,K] = computeMatrices(X,T,params,theReferenceElement)
%Matrices del sistema: masa Mlum Minv, rigidez K, amortigualmiento D
%e inextensibilidad C (3-tensor aplanado)
theta = params(1);
alfa = params(2);
beta = params(3);
%elemento
w = theReferenceElement.IPweights;
N=theReferenceElement.N;
Nxi=theReferenceElement.Nxi;
Neta=theReferenceElement.Neta;
nodos_elemento = size(theReferenceElement.nodesCoord,1);

%mallado
n_nodos = size(X,1);
n_elementos = size(T,1);

%creamos las matrices del sistema
M = spalloc(n_nodos,n_nodos,9*n_nodos);
L = spalloc(n_nodos,n_nodos,9*n_nodos);
Cu = spalloc(n_nodos^2,n_nodos,27*n_nodos);
Cv = spalloc(n_nodos^2,n_nodos,27*n_nodos);
Cuv = spalloc(n_nodos^2,n_nodos,27*n_nodos);

%Loop in elements
for i=1:n_elementos
    Te=T(i,:); %index of the nodes in the element
    Xe=X(Te,:); %coordinates of the nodes in the element  
    %elemental matrices
    Me = zeros(nodos_elemento,nodos_elemento); 
    Le = zeros(nodos_elemento,nodos_elemento); 
    Cu_e = zeros(nodos_elemento^2,nodos_elemento);
    Cv_e = zeros(nodos_elemento^2,nodos_elemento);
    Cuv_e = zeros(nodos_elemento^2,nodos_elemento);
    %Bucle en puntos de integracion
    for k=1:length(w)
        %calculo de las derivadas en el elemento
        N_k = N(k,:); Nxi_k = Nxi(k,:); Neta_k =  Neta(k,:);
        %vectores tangentes
        phi_xi = Nxi_k*Xe; phi_eta = Neta_k*Xe;
        dphi = [phi_xi; phi_eta]';
        %normal
        %nu = cross(phi_xi,phi_eta);
        %elemento de area
        E = (phi_xi)*(phi_xi)';
        G = (phi_eta)*(phi_eta)';
        F = (phi_xi)*(phi_eta)';
        dS = sqrt(abs(E*G - F^2))*w(k);
        %metrica
        m = [E F; 
             F G];
        coefs = m\[Nxi_k;Neta_k];
        Nxyz_k = dphi*coefs;
        Nx_k = Nxyz_k(1,:); Ny_k = Nxyz_k(2,:); Nz_k = Nxyz_k(3,:);
        %lo juntamos todo
        Me = Me + kron(N_k',N_k)*dS; %2-tensor masa
        Le = Le + (Nx_k'*Nx_k + Ny_k'*Ny_k + Nz_k'*Nz_k)*dS; %2-tensor laplaciano 
        %las tres condiciones para preservar la metrica
        Cu_e = Cu_e + kron(Nxi_k',kron(Nxi_k',N_k))*dS; %3-tensor 
        Cv_e = Cv_e + kron(Neta_k',kron(Neta_k',N_k))*dS; %3-tensor 
        Cuv_e = Cuv_e + 0.5*(kron(Nxi_k',kron(Neta_k',N_k))+kron(Neta_k',kron(Nxi_k',N_k)))*dS; %3-tensor 
    end    
    %assembly of elemental 2-tensors
    M(Te,Te) = M(Te,Te) + Me;
    L(Te,Te) = L(Te,Te) + Le;
    %3 tensor
    indices = combvectores(Te,Te); 
    ii = indices(1,:); jj = indices(2,:);
    TeTe = (jj - 1)*n_nodos + ii;
    Cu(TeTe,Te) = Cu(TeTe,Te) + Cu_e; 
    Cv(TeTe,Te) = Cv(TeTe,Te) + Cv_e; 
    Cuv(TeTe,Te) = Cuv(TeTe,Te) + Cuv_e; 
end

M_lumped = sum(M,1); 
Minv = sparse(diag(1./M_lumped)); 
C = struct('Cu',Cu*Minv,'Cv',Cv*Minv,'Cuv',Cuv*Minv);
%matrices
Mlum = sparse(blkdiag(diag(M_lumped),diag(M_lumped),diag(M_lumped)));
Kben = sparse(blkdiag(L'*Minv*L,L'*Minv*L,L'*Minv*L));
Minv = sparse(blkdiag(Minv,Minv,Minv)); 
D = alfa*Mlum + beta*Kben;
K = theta*Kben;



