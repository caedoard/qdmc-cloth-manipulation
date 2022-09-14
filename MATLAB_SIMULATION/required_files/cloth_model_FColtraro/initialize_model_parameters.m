function [model,Xin] = initialize_model_parameters(npx,npy)
% Taken from F.Coltraro

% DEFINITION OF SIMULATION PARAMETERS
model.dt = 0.01;

% This parameters are taken from F.Coltraro validation:
model.rho = 0.2; %densidad de masa
theta = 0; %resistencia  doblarse: rigidez;
alfa = 1.5;  %amortiguamiento oscilaciones lentas (~rozamiento con el aire)
beta = 0; %amortiguamiento oscilaciones rapidas

% model.rho = 0.10; %densidad de masa
% theta = 10; %resistencia  doblarse: rigidez;
% alfa = 0.15;  %amortiguamiento oscilaciones lentas (~rozamiento con el aire)
% beta = 4; %amortiguamiento oscilaciones rapidas

% MALLADO INICIAL
%se crea a partir de una parametrizacion de la forma (f1(x,y),f2(x,y),f3(x,y))
%donde (x,y) estan en el rectangulo [ax,bx]x[ay,by]
ax = -0.15; bx = 0.15; ay = -0.15; by = 0.15;
f1 = @(x,y) x; f2 = @(x,y) 0.75*sqrt(1 - x.^2) -0.1; f3 = @(x,y) y;
model.n_nodos = npx*npy; %numero de puntos
[Xin,model.T] = CreateMesh([ax,bx,ay,by],npx,npy,f1,f2,f3); 


% BOUNDARY
%deteccion de los puntos del borde para poder imponer isometria tambien alli
[Xb,Tb,model.nodos_borde] = GetBoundary(Xin,[ax,bx,ay,by],npx,npy,f1,f2,f3); 
model.n_nodos_bnd = size(Xb.Xb,1);
%interior nodes
nodes_int = setdiff(1:model.n_nodos,model.nodos_borde.nodes_bnd); 
%coordenadas de las esquinas
esquinas = [1, npx, npx*(npy-1)+1, npx*npy]; 

% NODOS A CONTROLAR
%damos las coordenadas de los nodos de los que vamos a fijar su trayectoria
controlados = [npx*(npy-1)+1, npx*npy]; %dos esquinas
model.coord_controlados = [controlados, controlados+model.n_nodos, controlados+(2*model.n_nodos)];
%matriz para imponer las condiciones en el contorno
model.A_b = spalloc(length(model.coord_controlados),3*model.n_nodos,length(model.coord_controlados));
model.A_b(:,model.coord_controlados) = eye(length(model.coord_controlados));

% ELEMENTOS DE REFERENCIA
%cuadrilateros bilineales
theReferenceElement = createReferenceElement();

% SYSTEM OF EQUATIONS
%calculo de las matrices del sistema: masa Mlum Minv, rigidez K, 
%amortigualmiento D e inextensibilidad C (3-tensor)
[model.C,model.n_conds,model.Mlum,model.Minv,model.Mcons,model.D,model.K] = computeMatrices(Xin,model.T,...
                               nodes_int,model.nodos_borde,esquinas,...
                               [theta,alfa,beta],...
                               theReferenceElement);
[model.Cphi0,~] = fun_C(Xin,model.C,model.Mcons,model.A_b,model.n_nodos,model.n_conds);
                           
% GRAVEDAD Y FUERZAS EXTERNAS
model.Fg = model.Mlum*reshape([zeros([model.n_nodos,1]),...
                   zeros([model.n_nodos,1]),...
          -9.8*model.rho*ones([model.n_nodos,1])],[3*model.n_nodos,1]);
      
model.A_b = full(model.A_b);

      
end