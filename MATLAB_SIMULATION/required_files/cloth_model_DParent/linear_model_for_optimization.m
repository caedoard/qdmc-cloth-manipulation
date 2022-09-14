function [xdot2] = linear_model_for_optimization(x,u,u_bef,param)
% It computes the equation x(k+1) = A·x(k)+B·u(k)+f_ext in the case the
% matrix A and f_ext are unkown (are a function of the optimization
% parameters)
%
% Author: David Parent, davidparentalonso@gmail.com
% Last review: 01/02/2021

k = param.stiffness;
b = param.damping;
ts = param.dt;
m = param.mass;
g = param.grav;
mat_x = param.mat_x;
mat_y = param.mat_y;
mat_z = param.mat_z;


% conn matrix indicates the connections between the nodes. They are
% numbered from bottom to top and from left to right
% THE FOLLOWING ADJACENCY MATRIX LOOKS WRONG
% conn   = -1.*[0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
%               1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0;
%               0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0;
%               0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0;
%               1 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0;
%               0 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0;
%               0 0 1 0 0 1 0 1 0 0 1 0 0 0 0 0;
%               0 0 0 1 0 0 1 0 0 0 0 1 0 0 0 0;
%               0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0;
%               0 0 0 0 0 1 0 0 1 0 1 0 0 1 0 0;
%               0 0 0 0 0 0 1 0 0 1 0 1 0 0 1 0;
%               0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 1;
%               0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%               0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0;
%               0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1;
%               0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
          
conn   = -1.*[0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
              1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0;
              0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0;
              0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0;
              1 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0;
              0 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0;
              0 0 1 0 0 1 0 1 0 0 1 0 0 0 0 0;
              0 0 0 1 0 0 1 0 0 0 0 1 0 0 0 0;
              0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0;
              0 0 0 0 0 1 0 0 1 0 1 0 0 1 0 0;
              0 0 0 0 0 0 1 0 0 1 0 1 0 0 1 0;
              0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 1;
              0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0;
              0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0;
              0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1;
              0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0];
         
% Action-reaction makes that the diagonal is the sum of the connections         
diago = diag(-1*sum(conn')); %transposada pq sumi les files i no columnes

% "Unitary force"
F = (diago+conn);

% Springs' force in each direction of the space. For each direction, an
% stiffness has been asigned.
Fp = [k(1).*F zeros(16) zeros(16);
      zeros(16) k(2).*F zeros(16);
      zeros(16) zeros(16) k(3).*F];
  
% Dampers' force in each direction of the space. For each direction, an
% stiffness has been asigned.  
Fv = [b(1).*F zeros(16) zeros(16);
      zeros(16) b(2).*F zeros(16);
      zeros(16) zeros(16) b(3).*F];
 
% The system equation is: 
% x(k+1) = x(k)+dt*v(k)
% v(k+1) = v(k)+dt*a(k), where a(k) is the acceleration of the node,
% that can be computed as:
% a(k) = (1/m)*F_total = (1/m)*( K*x(k) + b*v(k) )
% With the computed matrices, this can be written as:
% [pos_next] =         I·[pos] +        dt·I·[vel]  
% [vel_next] = (dt/m)·Fp·[pos] + (I+dt/m·Fv)·[vel]  
% that is computed as:
A = [eye(48)      ts.*eye(48) ;
    (ts/m).*Fp    eye(48)+(ts/m).*Fv];
 
% The upper corners are set to 0 because they are fixed
A(13,49:end) = zeros(1,48);
A(16,49:end) = zeros(1,48);
A(29,49:end) = zeros(1,48);
A(32,49:end) = zeros(1,48);
A(45,49:end) = zeros(1,48);
A(48,49:end) = zeros(1,48);
 
A(13+48,:) = zeros(1,96);
A(16+48,:) = zeros(1,96);
A(29+48,:) = zeros(1,96);
A(32+48,:) = zeros(1,96);
A(45+48,:) = zeros(1,96);
A(48+48,:) = zeros(1,96);

% Control matrix has only 1 in the position of the upper corners
B = [zeros(12,6);
    [1 0 0 0 0 0];
    zeros(2,6);
    [0 1 0 0 0 0];
    zeros(12,6);
    [0 0 1 0 0 0];
    zeros(2,6);
    [0 0 0 1 0 0];
    zeros(12,6);
    [0 0 0 0 1 0];
    zeros(2,6);
    [0 0 0 0 0 1];
    zeros(48,6)];

% Initial lengths matrices
mx = sum(mat_x')'; 
my = sum(mat_y')';
mz = sum(mat_z')';
grav = [zeros(48,1); zeros(32,1); -g*ones(16,1)];
l0_molles = [zeros(48,1); -k(1).*mx./m; -k(2).*my./m; -k(3).*mz./m];

f_ext = grav + l0_molles;
f_ext([13+48 16+48 29+48 32+48 45+48 48+48]) = zeros(6,1);

xdot2 = A*x + B*(u-u_bef) + ts.*f_ext;
end