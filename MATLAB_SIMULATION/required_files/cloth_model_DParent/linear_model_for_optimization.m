function [xdot2] = linear_model_for_optimization(x,u,u_bef,param)
    % Creates the matrices to express the system in linear state space:
    % x(k+1) = A*x(k) + B*u(k) + f_ext
    %
    % - Original Author: David Parent, davidparentalonso@gmail.com
    % - Modified by: Adrià Luque, adria.alados@gmail.com
    % - Last review: September 2021

    k = param.stiffness;
    b = param.damping;
    Ts = param.dt;
    m = param.mass;
    g = param.grav;

    nx = param.row;
    ny = param.col;
    blksz = 3*nx*ny;
    cctrl = param.coord_ctrl;

    % Connectivity matrix (connections between nodes)
    % Numbered from bottom to top and from left to right
    conn = zeros(nx*ny);
    for i=0:nx*ny-1 
        if(mod(i,ny) ~= 0), conn(i+1, i-1+1) = -1; end
        if(mod(i,ny) ~= ny-1), conn(i+1, i+1+1) = -1; end
        if(floor(i/ny) ~= 0), conn(i+1, i-ny+1) = -1; end
        if(floor(i/ny) ~= nx-1), conn(i+1, i+ny+1) = -1; end
    end
    % Remove connections on controlled nodes
    conn(cctrl(1),:) = 0;
    conn(cctrl(2),:) = 0;


    % Action-reaction makes that the diagonal is the sum of the connections
    diago = -diag(sum(conn,2)); % row sum

    % "Unitary force"
    F = diago + conn;

    % Springs' force in each direction of the space. For each direction, an
    % stiffness has been asigned.
    Fp = [k(1)*F zeros(nx*ny) zeros(nx*ny);
          zeros(nx*ny) k(2)*F zeros(nx*ny);
          zeros(nx*ny) zeros(nx*ny) k(3)*F];

    % Dampers' force in each direction of the space. For each direction, an
    % stiffness has been asigned.  
    Fv = [b(1)*F zeros(nx*ny) zeros(nx*ny);
          zeros(nx*ny) b(2)*F zeros(nx*ny);
          zeros(nx*ny) zeros(nx*ny) b(3)*F];

    % The system equation is: 
    % x(k+1) = x(k) + dt*v(k)
    % v(k+1) = v(k) + dt*a(k)
    % Where: a(k) = (1/m)*F_total = (1/m)*(K*x(k) + b*v(k))
    % Thus:
    % [pos_next] =         I*[pos] +        dt*I*[vel]  
    % [vel_next] = (dt/m)*Fp*[pos] + (I+dt/m·Fv)*[vel]  
    A = [eye(blksz)      Ts*eye(blksz);
         (Ts/m)*Fp    eye(blksz)+(Ts/m)*Fv];

    % The upper corners are set to 0 because they are fixed
    A(cctrl, blksz+1:end) = 0;
    A(cctrl+blksz, :) = 0;

    % Control matrix has only 1 in the position of the upper corners
    B = zeros(2*blksz,6);
    for i=1:length(cctrl)
        B(cctrl(i), i) = 1;
    end

    % Matrices of initial neighbor distance (pull)
    % Row=node -> If neighbors at same distance in one direction: equilibrium
    %          -> If edge node: pulled by only one neighbor in a direction
    %          -> With z_sum, center nodes can have pull in z
    mx = sum(param.mat_x, 2);
    my = sum(param.mat_y, 2);
    mz = sum(param.mat_z, 2);

    grav = zeros(2*blksz, 1);
    grav(end-nx*ny+1:end) = -g;

    a0_molles = [zeros(blksz,1);
                 -k(1)/m * mx; -k(2)/m * my; -k(3)/m * mz];

    f_ext = grav + a0_molles;
    f_ext(cctrl+blksz) = 0;
    xdot2 = A*x + B*(u-u_bef) + Ts.*f_ext;
end