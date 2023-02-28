function [mdl, X0] = initialize_nl_model(l_side,n_side,cpt,angle,dt)

%% Initial parameters
mdl = struct();
mdl.dt = dt;

% Cloth parameters
mdl.rho = 0.2;   % Density (0.1 silk, 0.2 shirt, 0.08 paper)
delta = mdl.rho; % Gravitational Mass (for aerodynamics)
theta = 0;       % Resistencia  doblarse: rigidez;
alfa = 1.5;      % Amortiguamiento oscilaciones lentas 
beta = 0;        % Amortiguamiento oscilaciones rapidas


%% Initial mesh
npx = n_side(1);
if length(n_side) == 1
    npy = npx;
else
    npy = n_side(2);
end
mdl.n_nodos = npx*npy; % Number of nodes/points
ax = -l_side/2; bx = l_side/2;
ay = -l_side/2; by = l_side/2; 
f1 = @(x,z) cpt(1) + x*cos(angle);
f2 = @(x,z) cpt(2) + x*sin(angle);
f3 = @(x,z) cpt(3) + z;
[X0, mdl.T] = CreateMesh([ax,bx,ay,by],npx,npy,f1,f2,f3);


%% Boundary
[Xtri,Ttri] = triangulateQuadMesh(X0, mdl.T);
TR = triangulation(Ttri,Xtri); 

% Boundary edges
Eb = freeBoundary(TR);
n_aristas = size(Eb,1);

% Boundary nodes
mdl.nodos_borde = Eb(:,1);
mdl.n_nodos_bnd = size(mdl.nodos_borde,1);

% Interior nodes
nodes_int = setdiff(1:mdl.n_nodos, mdl.nodos_borde);

% Matrices para imponer isometria en el borde:
%  imponemos que se preserve 
%  la longitud de las aristas del borde
I = [1:n_aristas,1:n_aristas,1:n_aristas,...
     1:n_aristas,1:n_aristas,1:n_aristas];
J = [Eb(:,2)', Eb(:,2)' + mdl.n_nodos, Eb(:,2)' + 2*mdl.n_nodos,...
     Eb(:,1)', Eb(:,1)' + mdl.n_nodos, Eb(:,1)' + 2*mdl.n_nodos];
mdl.bnd_aux = struct('Eb',Eb,'I',I,'J',J);

esquinas = [1, npx, npx*(npy-1)+1, npx*npy]; 


%% Corner nodes and coordinates
nodes_ctrl = [npx*(npy-1)+1, npx*npy]; % Two corners
mdl.coord_ctrl = [nodes_ctrl, nodes_ctrl+mdl.n_nodos, nodes_ctrl+(2*mdl.n_nodos)];
mdl.coord_lc = [1 npx 1+npx*npy npx*npy+npx 2*npx*npy+1 2*npx*npy+npx];

% Matrix to impose boundary conditions
mdl.A_b = spalloc(length(mdl.coord_ctrl),3*mdl.n_nodos,length(mdl.coord_ctrl));
mdl.A_b(:,mdl.coord_ctrl) = eye(length(mdl.coord_ctrl));


%% Reference elements
% Bilineal quad
theReferenceElement = createReferenceElement();


%% System of Equations

% System matrices: mass Mlum Minv, stiffness K, 
% Damping D, Inextensibility C (3-tensor)
[mdl.C,mdl.Mlum,mdl.Minv,mdl.D,mdl.K] = computeMatrices(X0,mdl.T,...
                                        [theta,alfa,beta],...
                                        theReferenceElement);

% Integrate foces implicitly
E = mdl.rho*mdl.Mlum + mdl.dt*mdl.D + (mdl.dt^2)*mdl.K; 
mdl.dE = decomposition(E);

% Select and process constraints
[mdl.C,mdl.n_conds] = processConstraints(mdl.C,mdl.n_nodos, ...
                                         nodes_int,esquinas);

% Initial values of the constraints                        
[mdl.Cphi0,~,~] = fun_C(X0(:), mdl.C, mdl.A_b, mdl.bnd_aux, ...
                        mdl.n_nodos, mdl.n_conds);
mdl.den = mdl.Cphi0;
mdl.den(abs(mdl.den) < 10^-6) = Inf;


%% Gravity and External Forces
mdl.Fg = sparse(mdl.Mlum*reshape([zeros([mdl.n_nodos,1]), ...
                zeros([mdl.n_nodos,1]), ...
                -9.8*delta*ones([mdl.n_nodos,1])],[3*mdl.n_nodos,1]));


end