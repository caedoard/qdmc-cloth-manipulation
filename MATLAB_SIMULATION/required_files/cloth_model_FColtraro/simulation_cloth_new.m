close all; clc; clear;

%% Input parameters
% Time
tf = 500;
dt = 0.02;

% Size
l_side = 0.3;
n_side = 10;
cpt = [0,0.5,l_side/2];


%% Initialize model
[SOM, X0] = initialize_nl_model(l_side,n_side,cpt,0,dt);
%posaux = X0;
%X0(:,1) = posaux(:,2);
%X0(:,2) = -posaux(:,1);


%% Trajectory of the Controlled nodes
times = dt*(0:tf);
nodes_ctrl = SOM.coord_ctrl(1:2);
u = zeros(3*length(nodes_ctrl),tf+1);
u(1:3,1) = X0(nodes_ctrl(1),:);
u(4:6,1) = X0(nodes_ctrl(2),:);
for tt=2:(tf+1)
    if times(tt) >= 2 && times(tt) <= 5
        v = [0;1;0;0;1;0];
        %v = [1;0;0;1;0;0];
        freq = 0.5;
        u(:,tt) = u(:,tt-1) + 0.015*cos(2*pi*freq*(times(tt)-2))*v;
    else
        u(:,tt) = u(:,tt-1);
    end 
end

% Reorder u to get [x1 x2 y1 y2 z1 z2]
n_cont = length(nodes_ctrl);
ind_x = 1:3:3*n_cont;
ind_y = ind_x + 1;
ind_z = ind_y + 1; 
u_old = u; % Copy

u((1:n_cont) + 0*n_cont,:) = u_old(ind_x,:); 
u((1:n_cont) + 1*n_cont,:) = u_old(ind_y,:); 
u((1:n_cont) + 2*n_cont,:) = u_old(ind_z,:); 


%% Simulation loop
% Initial state
phi0 = sparse(X0(:)); 
dphi0 = sparse(zeros([3*SOM.n_nodos,1]));

% Storage
phiPositions = cell([1,tf+1]);
phiPositions{1} = reshape(phi0,[SOM.n_nodos,3]);
phiVelocities = cell([1,tf+1]);
phiVelocities{1} = reshape(dphi0,[SOM.n_nodos,3]);

XX = zeros([SOM.n_nodos,tf+1]);
YY = zeros([SOM.n_nodos,tf+1]);
ZZ = zeros([SOM.n_nodos,tf+1]);
XX(:,1) = X0(:,1);
YY(:,1) = X0(:,2);
ZZ(:,1) = X0(:,3);

tStart = tic;
for tt=1:tf
    % Display time
    disp(tt)
    
    % Control action
    u_t = u(:,tt+1);
    
    % Solve ODE
    [phi,dphi] = simulate_cloth_step([phi0;dphi0], u_t, SOM);
    
    % Save states
    phiPositions{tt+1} = reshape(phi,[SOM.n_nodos,3]);
    phiVelocities{tt+1} = reshape(dphi,[SOM.n_nodos,3]); 
    
    XX(:,tt+1) = phiPositions{tt+1}(:,1); 
    YY(:,tt+1) = phiPositions{tt+1}(:,2); 
    ZZ(:,tt+1) = phiPositions{tt+1}(:,3);
    
    % Update
    phi0 = phi;
    dphi0 = dphi;   
end

total = toc(tStart);
disp('Tiempo total:')
disp(total)


%% Movie
peli = ''; vel_vid = 1; 
hacerPeli(phiPositions, SOM.T, [], peli, vel_vid, dt)



