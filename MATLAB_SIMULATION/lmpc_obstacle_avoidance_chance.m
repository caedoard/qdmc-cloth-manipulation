% Solve the MPC problem where the COM is a linear cloth model and the SOM
% is the Coltraro's cloth model.
%
% Author: David Parent, davidparentalonso@gmail.com
% Last review: 07/02/2021

clear; close all; clc;


addpath('required_files/cloth_model_FColtraro')
addpath('required_files/cloth_model_DParent')
addpath('/Applications/casadi-matlabR2015a-v3.5.5')
import casadi.*


%% Generate trajectory and place obstacle, 4 phases of motion
error_seeds = [];
chance = true;
error_trajs = [];
runtime_trajs = [];
n_y = 6;
store_state = [];
store_u = [];
t1 = 120;
t2 = 120;
t3 = 120; 
t4 = 120;
t_step = 0.0025;
phi_l_Traj = [(-.30: t_step: 0)', (0: t_step: .30)', zeros(t1 + 1, 1)];
phi_l_Traj = [phi_l_Traj;
              (0: t_step: 0.30)', .30*ones(t2 + 1, 1), zeros(t2 + 1, 1)];
phi_l_Traj = [phi_l_Traj;
              (0.30: t_step: 0.60)', .30*ones(t3 + 1, 1), zeros(t3 + 1, 1)];
phi_l_Traj = [phi_l_Traj;
              (0.60: t_step: 0.90)', (0.30:-t_step:0)', zeros(t4 + 1, 1)];

phi_r_Traj = [(0: t_step: 0.30)', (0: t_step: .30)', zeros(t1 + 1, 1)];
phi_r_Traj = [phi_r_Traj;
              (0.30: t_step: 0.60)', .30*ones(t2 + 1, 1), zeros(t2 + 1, 1)];
phi_r_Traj = [phi_r_Traj;
              (0.60: t_step: 0.90)', .30*ones(t3 + 1, 1), zeros(t3 + 1, 1)];
phi_r_Traj = [phi_r_Traj;
              (0.90: t_step: 1.20)', (0.30:-t_step:0)', zeros(t4 + 1, 1)];
plot3(phi_l_Traj(:, 1), phi_l_Traj(:, 2), phi_l_Traj(:, 3), 'LineWidth', 2);
hold on
plot3(phi_r_Traj(:, 1), phi_r_Traj(:,2), phi_r_Traj(:, 3),'LineWidth', 2);
xc=0.45; yc=0.15; zc=-0.1;    % coordinated of the center
L=0.20;                 % cube size (length of an edge)
alpha=0.4;           % transparency (max=1=opaque)

X = [0 0 0 0 0 1; 1 0 1 1 1 1; 1 0 1 1 1 1; 0 0 0 0 0 1];
Y = [0 0 0 0 1 0; 0 1 0 0 1 1; 0 1 1 1 1 1; 0 0 1 1 1 0];
Z = [0 0 1 0 0 0; 0 0 1 0 0 0; 1 1 1 0 1 1; 1 1 1 0 1 1];

C='blue';                  % unicolor

X = L*(X-0.5) + xc;
Y = L*(Y-0.5) + yc;
Z = L*(Z-0.5) + zc; 

fill3(X,Y,Z,C,'FaceAlpha',alpha);  
line([0, 0.35], [0.3, 0.25], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'black');

line([0.55, 0.9], [0.25, 0.3], 'LineStyle', '--', 'LineWidth', 2, 'COlor', 'black');
set(gca, 'TickLabelInterpreter', 'latex');

legend('$\mathbf r_{\textit{left}}$', '$\mathbf r_{\textit{right}}$', 'Interpreter', 'latex',...
'Location', 'north outside', 'NumColumns', 2);
axis equal
xlabel('$x$ [m]', 'Interpreter', 'latex')
        ylabel('$y$ [m]', 'Interpreter', 'latex')
        grid on
hold off

x0_l = phi_l_Traj(1, :);
x0_r = phi_r_Traj(1, :);
lCloth = norm(x0_l - x0_r);
cCloth = (x0_l + x0_r) ./ 2 + [0, 0, lCloth / 2];
dphi_corners = x0_r - x0_l;
aCloth = atan2(dphi_corners(2), dphi_corners(1));
disturbance = true;
% Define COM parameters(computation oriented model)
tot_runtimes = [];
tot_errors = [];
for seed = (1:10)
    rng(seed);
    n_nodes_edge = 4;
    alpha_err = 0.1;
    N = 25; % Number of control intervals
    Ms = 25;
    runtimes=[];
    errors=[];
    tottimes = [];
    for M = Ms
        store_state = [];
        store_u = [];
        COM = struct;
        COM.row = n_nodes_edge; COM.col = n_nodes_edge;
        COM.mass = 0.1;
        COM.grav = 9.8;
        COM.dt = 0.02;
        nr = COM.row;
        nc = COM.col;
        % parameters from calibration
        COM.stiffness = [-41.7998405537012  -39.1879707090799 -91.4967890838762]; 
        COM.damping = [-1.72864429464933   -1.33085664353071   -2.2554123678956]; 
        COM.z_sum = 0.0709542759848484 + 0.007;
        % Define the SOM
        nx = 7; ny = nx;
        coord_l = [1 nc 1+nr*nc nr*nc+nc 2*nr*nc+1 2*nr*nc+nc]; %coordenates of the linear model (of the lower corners)
        coord_nl = [1 nx 1+nx*ny nx*ny+nx 2*nx*ny+1 2*nx*ny+nx]; %coordenates of the non-linear model (of the lower corners)
        COM_nd_ctrl = [nr*(nc-1)+1, nr*nc];
        COM.coord_ctrl = [COM_nd_ctrl, COM_nd_ctrl+nr*nc, COM_nd_ctrl+2*nr*nc];  
        Ts = COM.dt;
        [SOM, pos] = initialize_nl_model(lCloth,nx,cCloth,aCloth, Ts);

        % Definite initial position of the nodes (defined here because we need it
        % to compute ext_force
        x_ini_SOM = [reshape(pos,[3*nx*ny 1]);zeros(3*nx*ny,1)]; %initial velocity=0
        u_SOM = x_ini_SOM(SOM.coord_ctrl);
        cloth_x = u_SOM([2 4 6]) - u_SOM([1 3 5]);
        cloth_y = [-cloth_x(2) cloth_x(1) 0]';
        cloth_z = cross(cloth_x,cloth_y);

        cloth_x = cloth_x/norm(cloth_x);
        cloth_y = cloth_y/norm(cloth_y);
        cloth_z = cloth_z/norm(cloth_z);
        Rcloth = [cloth_x cloth_y cloth_z];
        % Simulate some SOM steps to stabilize the NL model
%         warning('off','MATLAB:nearlySingularMatrix');
%         lastwarn('','');
%         [p_ini_SOM, ~] = simulate_cloth_step(x_ini_SOM,u_SOM,SOM);
%         [~, warnID] = lastwarn;
%         while strcmp(warnID, 'MATLAB:nearlySingularMatrix')
%             lastwarn('','');
%             x_ini_SOM = [p_ini_SOM; zeros(3*nx*ny,1)];
%             [p_ini_SOM, ~] = simulate_cloth_step(x_ini_SOM,u_SOM,SOM);
%             [~, warnID] = lastwarn;
%         end
%         warning('on','MATLAB:nearlySingularMatrix');
        [reduced_pos,~] = take_reduced_mesh(x_ini_SOM(1:3*nx*ny), ...
                                        x_ini_SOM(3*nx*ny+1:6*nx*ny), ...
                                        nx, COM.row);
        x_ini_COM = [reduced_pos; zeros(3*COM.row*COM.col,1)];
        % Rotate initial COM position to XZ plane
        RCloth_ini = [cos(aCloth) -sin(aCloth) 0; sin(aCloth) cos(aCloth) 0; 0 0 1];
        posCOM = reshape(x_ini_COM(1:3*nr*nc), [nr*nc,3]);
        posCOM_XZ = (RCloth_ini^-1 * posCOM')';

        % Initial position of the nodes
        COM.nodeInitial = lift_z(posCOM_XZ, COM);
        % Find initial spring length in each direction x,y,z
        [COM.mat_x, COM.mat_y, COM.mat_z] = compute_l0_linear(COM,0);

        % Find linear matrices
        [A, B, ext_force] = create_model_linear_matrices(COM);
        B_delta = zeros(2 * COM.row * COM.col * 3, 3);
        B_delta(COM.row * COM.col * 3 + 1:COM.row * COM.col * 4, 1) = COM.dt * 1/COM.mass * ones(COM.row * COM.col, 1);
        B_delta(COM.row * COM.col * 4 + 1:COM.row * COM.col * 5, 2) = COM.dt * 1/COM.mass * ones(COM.row * COM.col, 1);
        B_delta(COM.row * COM.col * 5 + 1:COM.row * COM.col * 6, 3) = COM.dt * 1/COM.mass * ones(COM.row * COM.col, 1);
        B_delta([COM.row * COM.col - COM.row + 1 + COM.row * COM.col * 3 COM.row * COM.col + COM.row * COM.col * 3], 1) = zeros(2,1);
        B_delta([2 * COM.row * COM.col - COM.row + 1 + COM.row * COM.col * 3 COM.row * COM.col * 2 + COM.row * COM.col * 3], 2) = zeros(2,1);
        B_delta([3 * COM.row * COM.col - COM.row + 1 + COM.row * COM.col * 3 COM.row * COM.col * 3 + COM.row * COM.col * 3], 3) = zeros(2,1);
        B_delta = sparse(B_delta);
        %% Start casADi optimization problem
        % Declare model variables
        n_d = 3;

        nr = COM.row;
        nc = COM.col;
        n_states = nr*nc*3*2; %should be 96 in a 4x4 model
        ubound = 50e-3; %0.8e-3, 50e-3
        P = SX.sym('P',n_states + n_d + n_y + 1,N*n_y+3); %Parameters in the optimization problem. It includes initial state, input and reference

        % Start with an empty NLP
        w=[]; %variables to optimize
        lbw = []; %lower bounds
        ubw = []; %upper bounds
        obj = 0; %objective function
        g= []; %constraints
        lbg = [];
        ubg = [];

        % Initial condition (initial state)
        Xk = P(1:n_states,1);
        U0 = P(1:6, 2);
        Q = 0.05;
        R = 1; %10, 20
        take_u=[];
        i=1;
        U = SX.sym('U',6, N);
        delta_u = [U(:, 1) - U0, diff(U, 1, 2)];
        w = U(:);
        Cc = P(n_states + n_d + 1: end - 1, 4:end)';
        bc = P(end, 4:end)';
        for k = 1:N
            % New NLP variable for the control
            %The bounds of the control input are very relevant!! If too big, the COM
            %allows it but the SOM not (the COM is a spring that from a time step
            %to the next one can be mathematically infinetely stretched). The
            %bounds should be approximately the maximum displacement expected.
            lbw = [lbw; -ubound*ones(6,1)]; 
            ubw = [ubw;  ubound*ones(6,1)];

            % Integrate till the end of the interval
            if k <= M
                inputs = [U(:, k);
                          P(n_states + 1:n_states + n_d, k)];
            else
                inputs = [zeros(6, 1);
                          P(n_states + 1:n_states + n_d, k)];
            end
            % Update objective function
            Xk_next = A*Xk + B*inputs(1:6) + (COM.dt).*ext_force + B_delta * inputs(6 + 1:end);
            obj = obj + (Xk_next(coord_l)- P(1:6,k+1+2))'*Q*(Xk_next(coord_l)- P(1:6,k+1+2));
            obj = obj + R*delta_u(:, k)'*delta_u(:, k);
            x_ctrl = Xk_next(COM.coord_ctrl);
            node_l = COM.coord_ctrl([1, 3, 5]);
            node_r = COM.coord_ctrl([2, 4, 6]);
            g = [g; sum((x_ctrl([2,4,6]) - x_ctrl([1,3,5])).^2) - lCloth^2 ];
            gbound = 0;
            lbg = [lbg; -gbound];
            ubg = [ubg;  gbound];
            Xk = Xk_next;
            Cc_curr = Cc((k - 1) * n_y + 1: k * n_y, :);
            bc_curr = bc((k - 1) * n_y + 1: k * n_y, :);
            output = Xk(coord_l);
            g = [g; Cc_curr * output - bc_curr];
            lbg = [lbg; -inf * ones(6, 1)];
            ubg = [ubg; zeros(6, 1)];
        end
        
        nlp_prob = struct('f', obj, 'x', w, 'g', g, 'p', P);
        opts = struct;
        opts.print_time = 0;
        opts.ipopt.print_level = 0; %0 to print the minimum, 3 to print the maximum
        opts.ipopt.warm_start_init_point = 'yes'; %warm start

        solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

        %----------------------------------------------
        % ALL OF THE ABOVE IS JUST A PROBLEM SET UP

        %%
        % THE SIMULATION LOOP SHOULD START FROM HERE
        %-------------------------------------------

        % Define trajectory to follow
        %scatter3(phired(1:nr*nc), phired(nr*nc + 1:nr*nc*2), phired(nr*nc*2 + 1:nr*nc*3));
        SOM_gravity = SOM.Fg;
        u_ini = x_ini_SOM(SOM.coord_ctrl);
        u_bef = u_ini;

        % Initialize things
        time_vector = (1:size(phi_l_Traj,1) + 1)' * COM.dt;
        store_state(:,end + 1) = x_ini_SOM;
        store_u(:,end + 1) = zeros(6,1);
        tT = 0;
        e = exp(1);
        wind_velocity = 2.4 * (e * ones(length(time_vector), 1)).^(-((time_vector - 5) ./ 2).^2);
        if disturbance
            wind_force = -0.5 * 1.225 * wind_velocity.^2 * 0.09; 
            noise_var = 0.005;
            sigma_epsilon_i = noise_var*eye(3);
            last_sigma = B_delta * sigma_epsilon_i * B_delta';
            sigmas = last_sigma;
            for i = (1:N - 1) % propagate the uncertainty through the autoregressive model
                last_sigma = A * last_sigma * A';
                addendum = B_delta * sigma_epsilon_i * B_delta';
                last_sigma = last_sigma + addendum;
                sigmas = [sigmas; last_sigma];
            end
        else
            wind_force = 0 * wind_velocity; % acts along the y axis
            sigmas = zeros(n_states*N, n_states);
            noise_var = 0;
        end
        plot(wind_force);
        %t0 = tic;
        delta_u_bef = zeros(6, 1);
        for i=2:size(phi_l_Traj,1)
            %SOM_gravity = SOM.Fg;
            if i>=size(phi_l_Traj,1)-N-1 %the last N timesteps, trajectory should remain constant
                Traj_l_Hp = repmat(phi_l_Traj(end,:),N+1, 1);
                Traj_r_Hp = repmat(phi_r_Traj(end,:),N+1, 1);
            else
                Traj_l_Hp = phi_l_Traj(i:i+N,:);
                Traj_r_Hp = phi_r_Traj(i:i+N,:);
            end
            pos_ini_COM = reshape(x_ini_COM(1:3*nr*nc),[nr*nc,3]);
            vel_ini_COM = reshape(x_ini_COM(3*nr*nc+1:6*nr*nc),[nr*nc,3]);

            pos_ini_COM_rot = (Rcloth^-1 * pos_ini_COM')';
            vel_ini_COM_rot = (Rcloth^-1 * vel_ini_COM')';
            x_ini_COM_rot = [reshape(pos_ini_COM_rot,[3*nr*nc,1]);
            reshape(vel_ini_COM_rot,[3*nr*nc,1])];
            % Rotate reference trajectory to cloth base
            Traj_l_Hp_rot = (Rcloth^-1 * Traj_l_Hp')';
            Traj_r_Hp_rot = (Rcloth^-1 * Traj_r_Hp')';
%             scatter3(x_ini_COM_rot(1:nr*nc), x_ini_COM_rot(nr*nc + 1:nr*nc*2), x_ini_COM_rot(nr*nc*2 + 1:nr*nc*3));
%             
%             hold on
%             %scatter3(phi_ini_SOM(1:nx*ny), phi_ini_SOM(nx*ny + 1:nx*ny*2), phi_ini_SOM(nx*ny*2 + 1:nx*ny*3));
%             scatter3(Traj_l_Hp_rot(1,1), Traj_l_Hp_rot(1,2), Traj_l_Hp_rot(1,3));
%             scatter3(Traj_r_Hp_rot(1,1), Traj_r_Hp_rot(1,2), Traj_r_Hp_rot(1,3));
%             hold off
%             ylim([0, 2.0]);
%             xlim([-1, 1]);
%             zlim([-0.5, 1.5]);

            % Define reference in the prediction horizon (moving window)
            reference = zeros(6*COM.row*COM.col,N+1+2);
            reference(:,1) = x_ini_COM_rot;
            reference(:,2) = [delta_u_bef; zeros(n_states - 6, 1)];
            reference(1,3:end) = Traj_l_Hp_rot(:,1)';
            reference(2,3:end) = Traj_r_Hp_rot(:,1)';
            reference(3,3:end) = Traj_l_Hp_rot(:,2)';
            reference(4,3:end) = Traj_r_Hp_rot(:,2)';
            reference(5,3:end) = Traj_l_Hp_rot(:,3)';
            reference(6,3:end) = Traj_r_Hp_rot(:,3)';
            sliced_wind_force = wind_force(i:min(i+N-1, end))';
            sliced_wind_force = [sliced_wind_force, zeros(1, max(0, N +2 +1- length(sliced_wind_force)))];
            predicted_disturbance = [-1/3 * sliced_wind_force;
                                     sliced_wind_force;
                                     zeros(1, N+2+1)];

            % Initial seed of the optimization (for "u")
            args.x0 = [reshape(diff(reference(1:6,3:N+2),1,2), 6*(N-1), 1); zeros(6, 1)];
            Cconstr = [];
            bconstr = [];
            for k = 0:N - 1
                if i < t1 - 1
                    Ccurr = [[eye(2, 2), zeros(2, 4)];
                               zeros(4, 6)];
                    bcurr = [repmat(0.35, 2, 1);
                               zeros(4, 1)];
                elseif i >= t1 - 1 && i <= t1 + t2 + t3 + 1
                    Ccurr = [zeros(2, 6);
                            [zeros(2, 2), -1*eye(2, 2), zeros(2, 2)];
                             zeros(2, 6)];
                    bcurr = [zeros(2, 1);
                               repmat(-0.25, 2, 1);
                               zeros(2, 1)];
                elseif i > t1 + t2 + t3 + 1
                    Ccurr = [[-1*eye(2, 2), zeros(2, 4)];
                               zeros(4, 6)];
                    bcurr = [repmat(-0.55, 2, 1);
                               zeros(4, 1)];
                else
                    Ccurr = zeros(6);
                    bcurr = zeros(6, 1);
                end
                curr_sigma = sigmas(k*n_states + 1: (k+1) * n_states, :);
                curr_sigma_outputs = curr_sigma(coord_l, coord_l);
                curr_sigma_outputs = Ccurr * curr_sigma_outputs * Ccurr';
                correction = norminv(1 - alpha_err/(2 * n_y * N)) * diag(curr_sigma_outputs).^0.5;
                bcurr = bcurr - correction;
                Cconstr = [Cconstr; Ccurr];
                bconstr = [bconstr; bcurr];
            end
            
            Cconstr = sparse([zeros(3, n_y); Cconstr]);
            bconstr = sparse([zeros(3, 1); bconstr]);
            reference = [reference, zeros(length(reference), N * n_y - N)];
            predicted_disturbance = [predicted_disturbance, zeros(n_d, N * n_y - N)];
            args.p = [reference;
                      predicted_disturbance;
                      Cconstr';
                      bconstr'];

            % Find the solution "sol"
            t0 = tic;
            sol = solver('x0', args.x0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg, 'p', args.p);
            u_rot =  reshape(full(sol.x),6,N);
            u_sliced = u_rot(:, 1);
            delta_u_bef = u_sliced; % Optimization variable is displacement delta_u; then we penalize the upper corner accel.
            u_rot_l = u_sliced([1 3 5], 1);
            u_rot_r = u_sliced([2 4 6], 1);
            u_base = Rcloth * [u_rot_l, u_rot_r];
            u = reshape(u_base',[6,1]);
            stats = solver.stats();
            if stats.success == 0
                "Optimization failed"
                u_SOM = u_SOM;
            else
                u_SOM = u + u_bef;
            end
            u_bef = u_SOM; %update

            tT=tT+toc(t0);
            SOM_perturbation = [1/3 * repmat(-sliced_wind_force(1) + normrnd(0, sqrt(noise_var)), nx*ny, 1);
                                   repmat(sliced_wind_force(1) - normrnd(0, sqrt(noise_var)), nx*ny, 1);
                                   zeros(nx * ny, 1)];
            %hold on
            SOM_perturbation = SOM_perturbation / COM.mass;
            SOM_perturbation = SOM.Mlum*SOM_perturbation*SOM.rho;
            SOM.Fg = SOM_gravity + SOM_perturbation; % Wind acting on the y axis.

            [next_pos_SOM, next_vel_SOM] = simulate_cloth_step(store_state(:,i-1),u_SOM,SOM);

            phi_ini_SOM = full(next_pos_SOM);
            dphi_ini_SOM = full(next_vel_SOM);    
            %t0 = tic;
            cloth_x = u_SOM([2 4 6]) - u_SOM([1 3 5]);
            cloth_y = [-cloth_x(2) cloth_x(1) 0]';
            cloth_z = cross(cloth_x,cloth_y);

            cloth_x = cloth_x/norm(cloth_x);
            cloth_y = cloth_y/norm(cloth_y);
            cloth_z = cloth_z/norm(cloth_z);
            Rcloth = [cloth_x cloth_y cloth_z];
                
            % Close the loop
            [phired,dphired] = take_reduced_mesh(phi_ini_SOM,dphi_ini_SOM, nx, COM.row);

            x_ini_COM = [phired;dphired];

%             scatter3(x_ini_COM(1:nr*nc), x_ini_COM(nr*nc + 1:nr*nc*2), x_ini_COM(nr*nc*2 + 1:nr*nc*3));
%             hold on
%             ylim([-0.4, 0.7]);
%             xlim([-0.4, 0.8]);
%             zlim([-0.2, 0.6]);
%             hold off
            % Store things
            store_state(:,end + 1) = [phi_ini_SOM;dphi_ini_SOM];
            store_u(:,end + 1) = u(1);    

            if(mod(i,100)==0)
                fprintf(['Iter: ', num2str(i), ...
                         ' \t Avg. time/iter: ', num2str(tT/i*1000), ' ms \n']);
            end
        end
        tT = tT+toc(t0);

        %% Errors                 
        error_l = store_state(coord_nl([1,3,5]),:)'-phi_l_Traj;
        error_r = store_state(coord_nl([2,4,6]),:)'-phi_r_Traj;
        error_p = [vecnorm(error_l,2,2), vecnorm(error_r,2,2)];

        eMAE = mean(abs([error_l error_r]));
        eRMSE = sqrt(mean([error_l error_r].^2));

        eRMSEp = mean([norm(eRMSE([1,3,5]),2) norm(eRMSE([2,4,6]),2)]);
        eRMSEm = mean(eRMSE,2);
        eMAEm  = mean(eMAE,2); % Old "avg_error"
        fprintf(['-----------------------------------------\n', ...
                 ' -- Total time: \t',num2str(tT),' s \n', ...
                 ' -- Avg. t/iter:\t',num2str(tT/size(phi_l_Traj,1)*1000),' ms \n', ...
                 ' -- Avg. error: \t',num2str(1000*eMAEm),' mm \n', ...
                 ' -- Norm RMSE:  \t',num2str(1000*eRMSEp),' mm\n']);


        %% PLOTS
            time = 0:0.01:size(store_state,2)*0.01-0.01;
% 
%             fig1 = figure(1);
%             fig1.Color = [1,1,1];
%             fig1.Units = 'normalized';
%             fig1.Position = [0.5 0.6 0.5 0.3];
% 
%             subplot(7,2,1:2:12);
%             pa1som=plot(time,store_state(coord_nl([1 3 5]),:)','linewidth',1.5);
%             hold on
%             pa1ref=plot(time,phi_l_Traj(:,:),'--k','linewidth',1.2);
%             hold off
%             title('\textbf{Lower Left Corner}', 'Interpreter', 'latex')
%             grid on
%             xlabel('Time [s]', 'Interpreter', 'latex')
%             ylabel('Position [m]', 'Interpreter', 'latex')
%             xlim([0 round(time(end))])
%             set(gca, 'TickLabelInterpreter', 'latex');
% 
%             subplot(7,2,2:2:12);
%             plot(time,store_state(coord_nl([2 4 6]),:)','linewidth',1.5)
%             hold on
%             plot(time,phi_r_Traj(:,:)','--k','linewidth',1.2)
%             hold off
%             title('\textbf{Lower Right Corner}', 'Interpreter', 'latex')
%             grid on
%             xlabel('Time [s]', 'Interpreter', 'latex')
%             ylabel('Position [m]', 'Interpreter', 'latex')
%             xlim([0 round(time(end))])
%             set(gca, 'TickLabelInterpreter', 'latex');
% 
%             Lgnd1 = legend([pa1som' pa1ref(1)], ...
%                            '$p_x\ (y_1,y_2)$','$p_y\ (y_3,y_4)$', '$p_z\ (y_5,y_6)$', '$r$', ...
%                            'Orientation','horizontal', 'Interpreter', 'latex');
%             Lgnd1.Position(1) = 0.5-Lgnd1.Position(3)/2;
%             Lgnd1.Position(2) = 0.02;
%             Lgnd1.Position(3:4) = Lgnd1.Position(3:4) + 0.01;
        runtimes(end+1) = tT/size(phi_l_Traj,1)*1000;
        errors(end+1) = 1000*eMAEm;
        tottimes(end+1) = tT;
    end
    tot_runtimes = [tot_runtimes; runtimes];
    tot_errors = [tot_errors; errors];
end

writematrix(tot_errors, "obstacle_lmpc_error_sim_traj" + ".csv");
writematrix(tot_runtimes, "obstacle_lmpc_time_sim_traj" +".csv");

%     writematrix(errors, "error_sim_traj" + traj_index +".csv");
%     writematrix(runtimes, "time_sim_traj" + traj_index +".csv");
%     writematrix(tottimes, "tot_time_sim_traj" + traj_index +".csv");
