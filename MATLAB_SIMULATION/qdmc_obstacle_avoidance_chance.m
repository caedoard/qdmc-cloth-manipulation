% Solve the QDMC problem where the COM is a linear cloth model and the SOM
% is Coltraro's cloth model.
%
% Author: Edoardo Caldarelli, ecaldarelli@iri.upc.edu
% September 2022

clear; close all; clc;

addpath('required_files/cloth_model_FColtraro')
addpath('required_files/cloth_model_DParent')
addpath('/Applications/casadi-matlabR2015a-v3.5.5')
import casadi.*

%% COM definition
error_seeds = [];
error_trajs = [];
runtime_trajs = [];

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
plot(phi_l_Traj(:, 1), phi_l_Traj(:, 2), 'LineWidth', 2);
hold on
plot(phi_r_Traj(:, 1), phi_r_Traj(:,2),'LineWidth', 2);
x = [0.35 0.55 0.55 0.35];
y = [0.25 0.25 0.05 0.05];
alpha = 0.4;
ylim([-0.1, 0.31]);
line([0.35, 0.55], [0.25, 0.25], 'LineWidth', 1.5, 'Color', 'blue');
line([0.55, 0.55], [0.25, 0.05], 'LineWidth', 1.5, 'Color', 'blue');
line([0.55, 0.35], [0.05, 0.05], 'LineWidth', 1.5, 'Color', 'blue');
line([0.35, 0.35], [0.05, 0.25], 'LineWidth', 1.5, 'Color', 'blue');
fill(x,y,'b', 'FaceAlpha', alpha)
line([0, 0.35], [0.3, 0.25], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'black');

line([0.55, 0.9], [0.25, 0.3], 'LineStyle', '--', 'LineWidth', 2, 'COlor', 'black');
set(gca, 'TickLabelInterpreter', 'latex');

legend('$\mathbf r_{\textit{left}}$', '$\mathbf r_{\textit{right}}$', 'Interpreter', 'latex',...
'Location', 'north outside', 'NumColumns', 2);
%axis equal
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
% Define COM parameters (computation oriented model)
tot_runtimes = [];
tot_errors = [];
for chance = [true, false]
    for alpha_err = [1.0, 0.5, 0.1]
        rng(1);
        n_nodes_edge = 4;
        N = 500; % finite-step-response horizon until steady state is reached.
        H = 25;
        Ms = 13; % Control horizon, set it to 5 to avoid the tuning of this parameter
        runtimes = [];
        errors = [];
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

            coord_l = [1 nc 1+nr*nc nr*nc+nc 2*nr*nc+1 2*nr*nc+nc]; %coordinates of the linear model (of the lower corners)
            % parameters from calibration    
            COM.stiffness = [-41.7998405537012  -39.1879707090799 -91.4967890838762]; 
            COM.damping = [-1.72864429464933   -1.33085664353071   -2.2554123678956]; 
            COM.z_sum = 0.0709542759848484 + 0.007;
            % Define the SOM
            nx = 7; ny = nx;
            coord_nl = [1 nx 1+nx*ny nx*ny+nx 2*nx*ny+1 2*nx*ny+nx]; %coordinates of the non-linear model (of the lower corners)
            COM_nd_ctrl = [nr*(nc-1)+1, nr*nc];
            COM.coord_ctrl = [COM_nd_ctrl, COM_nd_ctrl+nr*nc, COM_nd_ctrl+2*nr*nc];
            % Definite initial position of the nodes (defined here because we need it
            % to compute ext_force
            Ts = COM.dt;
            [SOM, pos] = initialize_nl_model(lCloth,nx,cCloth,aCloth, Ts);
            x_ini_SOM = [reshape(pos,[3*nx*ny 1]);zeros(3*nx*ny,1)]; %initial velocity=0

            u_SOM = x_ini_SOM(SOM.coord_ctrl);

            cloth_x = u_SOM([2 4 6]) - u_SOM([1 3 5]);
            cloth_y = [-cloth_x(2) cloth_x(1) 0]';
            cloth_z = cross(cloth_x,cloth_y);

            cloth_x = cloth_x/norm(cloth_x);
            cloth_y = cloth_y/norm(cloth_y);
            cloth_z = cloth_z/norm(cloth_z);
            Rcloth = [cloth_x cloth_y cloth_z];
            warning('off','MATLAB:nearlySingularMatrix');
            lastwarn('','');
            [p_ini_SOM, ~] = simulate_cloth_step(x_ini_SOM,u_SOM,SOM);
            [~, warnID] = lastwarn;
    %             while strcmp(warnID, 'MATLAB:nearlySingularMatrix')
    %                 lastwarn('','');
    %                 x_ini_SOM = [p_ini_SOM; zeros(3*nx*ny,1)];
    %                 [p_ini_SOM, ~] = simulate_cloth_step(x_ini_SOM,u_SOM,SOM);
    %                 [~, warnID] = lastwarn;
    %                 scatter3(x_ini_SOM(1: nx*ny, :), x_ini_SOM(nx*ny + 1: 2 * nx*ny, :), x_ini_SOM(2*nx*ny + 1: 3*nx*ny, :));
    % 
    %             end
            warning('on','MATLAB:nearlySingularMatrix');
            x_ini_SOM(3*nx*ny + 1: end, :) = zeros(3*nx*ny, 1);
            warning('on','MATLAB:nearlySingularMatrix');
            [reduced_pos,~] = take_reduced_mesh(x_ini_SOM(1:3*nx*ny), ...
                                            x_ini_SOM(3*nx*ny+1:6*nx*ny), ...
                                            nx, COM.row);
    %         scatter3(reduced_pos(1:nr*nc), reduced_pos(nr*nc + 1:nr*nc*2), reduced_pos(nr*nc*2 + 1:nr*nc*3));
    %         hold on
    %         ylim([0.5, 0.7]);
    %         xlim([-0.2, 0.4]);
    %         zlim([-0.2, 0.3]);
    %         hold off
            x_ini_COM = [reduced_pos; zeros(3*COM.row*COM.col,1)];
            pos_ini_COM = reshape(x_ini_COM(1:3*nr*nc),[nr*nc,3]);
            vel_ini_COM = reshape(x_ini_COM(3*nr*nc+1:6*nr*nc),[nr*nc,3]);

            pos_ini_COM_rot = (Rcloth^-1 * pos_ini_COM')';
            vel_ini_COM_rot = (Rcloth^-1 * vel_ini_COM')';
            x_ini_COM_rot = [reshape(pos_ini_COM_rot,[3*nr*nc,1]); reshape(vel_ini_COM_rot,[3*nr*nc,1])];
            % Initial position of the COM nodes
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

            % Build disturbance matrix (wind is 3-dimensional disturbance)
            B_delta = zeros(2 * COM.row * COM.col * 3, 3);
            B_delta(COM.row * COM.col * 3 + 1:COM.row * COM.col * 4, 1) = COM.dt * 1/COM.mass * ones(COM.row * COM.col, 1);
            B_delta(COM.row * COM.col * 4 + 1:COM.row * COM.col * 5, 2) = COM.dt * 1/COM.mass * ones(COM.row * COM.col, 1);
            B_delta(COM.row * COM.col * 5 + 1:COM.row * COM.col * 6, 3) = COM.dt * 1/COM.mass * ones(COM.row * COM.col, 1);
            B_delta([COM.row * COM.col - COM.row + 1 + COM.row * COM.col * 3 COM.row * COM.col + COM.row * COM.col * 3], 1) = zeros(2,1);
            B_delta([2 * COM.row * COM.col - COM.row + 1 + COM.row * COM.col * 3 COM.row * COM.col * 2 + COM.row * COM.col * 3], 2) = zeros(2,1);
            B_delta([3 * COM.row * COM.col - COM.row + 1 + COM.row * COM.col * 3 COM.row * COM.col * 3 + COM.row * COM.col * 3], 3) = zeros(2,1);
            B_delta = sparse(B_delta);
            %% Step response for FSR model

            % Compute the coefficients of the step response of MIMO system.
            % Note that the system considered is linear, thus it suffices to 
            % apply one step for each input-output pair.

            n_u = 6;
            n_y = 6;
            step_amplitude = 1e-5; % Offset w.r.t. initial position
            step_response = zeros(N, length(x_ini_COM));
            free_response = zeros(10 * N, length(x_ini_COM));

            MIMO_step_response = zeros(0, N, n_y); % MIMO_step_response(a, b, c) is the output response of channel c at timestep b
                                                   % given a step input in channel a

            % WE MUST START AT THE EQUILIBRIUM
            xnew = x_ini_COM_rot;
            %xnew = A * xnew + B * [0.05; 0; 0; -0.5; 0; 0] + COM.dt .* ext_force;
            for k = (1: 15 * N)
                xnew = A * xnew + COM.dt .* ext_force;
                free_response(k, :) = xnew';
            end
            x_eq = xnew;
            free_output = free_response(:, coord_l);
            % For each input 
            counter = 1;
            for h = (1:n_u)
                delta_u_step = zeros(N, n_u);
                delta_u_step(1, h) = step_amplitude;
                % Compute step response
                xnew = x_eq;
                for k = (1:N)
                    xnew = A * xnew + B * delta_u_step(k, :)' + COM.dt .* ext_force;
                    step_response(k, :) = xnew';
                end
                % Extract output
                output_response = step_response(:, coord_l);
                step_coefficients = (output_response - x_eq(coord_l)') ./ step_amplitude;
                MIMO_step_response(h, :, :)  = step_coefficients;
    %                         pp = subplot(3,3,h);
    %             format = '%i';
    %             val = sprintf(format, h);
    %             title_string = "\textbf{Response coefficients for step input in channel " + val + "}";
    %             title(title_string, 'Interpreter', 'latex')
    %     
    %             for k = (1:6)
    %                 grid on
    %                 hold on
    %                 plot((1:N) .* COM.dt, step_coefficients(:, k), 'linewidth',1.5);
    %             end
    %             hold off
    %             legend("$S_{" + h + ", 1}$","$S_{" + h + ", 2}$", "$S_{" + h + ", 3}$",...
    %                    "$S_{" + h + ", 4}$", "$S_{" + h + ", 5}$", "$S_{" + h + ", 6}$", ...
    %                    'Interpreter', 'latex', 'Orientation','horizontal', 'location', 'southoutside');
    %     
    %             xlabel('t [s]', 'Interpreter', 'latex');
    %     
    %             set(gca, 'TickLabelInterpreter', 'latex');
            end

            % Now, identify the transfer from each disturbance channel (3) to each
            % output channel
            %figure
            n_d = 3;
            step_amplitude = 1e0; % Offset w.r.t. initial position
            step_response = zeros(N, length(x_ini_COM_rot));
            free_response = zeros(10 * N, length(x_ini_COM_rot));

            MIMO_step_response_disturbance = zeros(0, N, n_y); % MIMO_step_response(a, b, c) is the output response of channel c at timestep b
                                                   % given a step input in channel a

            % WE MUST START AT THE EQUILIBRIUM
            xnew = x_ini_COM_rot;
            for k = (1: 10 * N)
                xnew = A * xnew + COM.dt .* ext_force;
                free_response(k, :) = xnew';
            end
            x_eq = xnew;
            free_output = free_response(:, coord_l);
            % For each disturbance channel 
            val=1;
            for h = (1:n_d)
                d_step = zeros(N, n_d);
                d_step(:, h) = step_amplitude;
                % Compute step response
                xnew = x_eq;
                for k = (1:N)
                    xnew = A * xnew + B_delta * d_step(k, :)' + COM.dt .* ext_force;
                    step_response(k, :) = xnew';
                end
                % Extract output
                output_response = step_response(:, coord_l);
                step_coefficients = (output_response - x_eq(coord_l)') ./ step_amplitude;
                MIMO_step_response_disturbance(h, :, :)  = step_coefficients;
                val = val + 1;
            end

            %% QDMC setup

            % Rearrange matrix S containing all step response coefficients from each
            % input to each output, in a lower-triangular fashion.

            S_tilde = zeros(n_y * H, M * n_u);

            for i = (1:n_u)
                for j = (1:n_y)
                    % Construct lower-triangular matrix S_ij by shifting the step
                    % response vector accoringly.
                    S_ij_vec = MIMO_step_response(i, 1:H, j);
                    S_ij_mat = zeros(H, M);
                    for h = (1:M)
                        S_ij_mat(h:end, h) = S_ij_vec(1:end + 1 - h);
                    end
            %         figure
            %         plot(S_ij_mat(:, 1))
                    S_tilde((j - 1) * H + 1:j * H, (i - 1) * M + 1: i * M) = S_ij_mat;
                end
            end

            S_tilde_steady_state = zeros(n_y * N, M * n_u);

            for i = (1:n_u)
                for j = (1:n_y)
                    % Construct lower-triangular matrix S_ij by shifting the step
                    % response vector accoringly.
                    S_ij_vec = MIMO_step_response(i, 1:N, j);
                    S_ij_mat = zeros(N, M);
                    for h = (1:M)
                        S_ij_mat(h:end, h) = S_ij_vec(1:end + 1 - h);
                    end
                    S_tilde_steady_state((j - 1) * N + 1:j * N, (i - 1) * M + 1: i * M) = S_ij_mat;
                end
            end

            % Now build the response matrix for the disturbance

            D_tilde = zeros(n_y * H, H * n_d);

            for i = (1:n_d)
                for j = (1:n_y)
                    % Construct lower-triangular matrix S_ij by shifting the step
                    % response vector accoringly.
                    D_ij_vec = MIMO_step_response_disturbance(i, 1:H, j);
                    D_ij_mat = zeros(H, H);
                    for h = (1:H)
                        D_ij_mat(h:end, h) = D_ij_vec(1:end + 1 - h);
                    end
                    D_tilde((j - 1) * H + 1:j * H, (i - 1) * H + 1: i * H) = D_ij_mat;
                end
            end

            D_tilde_steady_state = zeros(n_y * N, H * n_d);

            for i = (1:n_d)
                for j = (1:n_y)
                    % Construct lower-triangular matrix S_ij by shifting the step
                    % response vector accoringly.
                    D_ij_vec = MIMO_step_response_disturbance(i, 1:N, j);
                    D_ij_mat = zeros(N, H);
                    for h = (1:H)
                        D_ij_mat(h:end, h) = D_ij_vec(1:end + 1 - h);
                    end
                    D_tilde_steady_state((j - 1) * N + 1:j * N, (i - 1) * H + 1: i * H) = D_ij_mat;
                end
            end

            % Compute initial free step response. Define
            % shifting matrix as well (needed to update the free step response at each 
            % time-step).
            Gamma = [];
            for i = (1:n_y)
                curr_shifting_matrix = eye(N - 1, N - 1);
                curr_shifting_matrix = [zeros(N - 1, 1), curr_shifting_matrix];
                curr_shifting_matrix = [curr_shifting_matrix;
                                        curr_shifting_matrix(end, :)];
                left_zero_padding = zeros(N, (i-1) * N);
                right_zero_padding = zeros(N, (n_y - i) * N);
                padded_shifting_matrix = [left_zero_padding, curr_shifting_matrix, right_zero_padding];
                Gamma = [Gamma;
                    padded_shifting_matrix];
            end

            % Define the shifting matrix needed for predicting y_{k} = Psi * f_{k} +
            % .... y_{k} = [y_{k+1|k}, y_{k+2|k}, ...] depends on [f_{k+1|k},
            % f_{k+2|k}, ...], not on [f_{k+1|k+1}, f_{k+2|k+1}, ...]

            Psi = [];
            for i = (1:n_y)
                curr_shifting_matrix = eye(H, H);
                curr_shifting_matrix = [zeros(H, 1), curr_shifting_matrix];
                curr_shifting_matrix = [curr_shifting_matrix, zeros(H, N - H - 1)];
                left_zero_padding = zeros(H, (i-1) * N);
                right_zero_padding = zeros(H, (n_y - i) * N);
                padded_shifting_matrix = [left_zero_padding, curr_shifting_matrix, right_zero_padding];
                Psi = [Psi;
                    padded_shifting_matrix];
            end
            S_tilde(S_tilde<1e-8) = 0;
            S_tilde_steady_state(S_tilde_steady_state<1e-8) = 0;
            D_tilde(D_tilde<1e-12) = 0;
            D_tilde_steady_state(D_tilde_steady_state<1e-8) = 0;
            Psi = sparse(Psi);
            Gamma = sparse(Gamma);
            S_tilde = sparse(S_tilde);
            S_tilde_steady_state = sparse(S_tilde_steady_state);
            D_tilde = sparse(D_tilde);
            D_tilde_steady_state = sparse(D_tilde_steady_state);
            % Build optimization risk (casADi symbolic objective function).
            % Declare model variables
            P = SX.sym('P', H * n_y, 4 + 6 * H + 1); % Parameters. First col is f_tilde_k, second 
                                          % col is reference signal, third col is
                                          % predicted disturbance, fourth col
                                          % is previous delta ; we also
                                          % pass the matrices for the
                                          % chance constraints!
            ubound = 50e-3; %0.8e-3, 50e-3
            % Tune R and Q
            % error
            R = 1; 
            Q = 0.01;
            % Start with an empty NLP
            w=[]; %variables to optimize
            lbw = []; %lower bounds
            ubw = []; %upper bounds
            obj = 0; %objective function
            g = [];
            lbg = [];
            ubg = [];
            writematrix(full(Psi), "Psi.csv");
            writematrix(full(D_tilde), "D_tilde.csv");
            writematrix(full(D_tilde_steady_state), "D_tilde_ss.csv");
            writematrix(full(S_tilde), "S_tilde.csv");
            writematrix(full(S_tilde_steady_state), "S_tilde_ss.csv");
            writematrix(full(Gamma), "Gamma.csv");
            writematrix(full(B_delta), "B_delta.csv");


            %% Symbolic risk assembly at tstep k *of the MPC loop*
            % New NLP variable for the control
            U = SX.sym('U', n_u , M);
            U0 = P(1:6, 4);
            Cc = P(:, 4 + 1: end - 1);
            bc = P(:, end);
            %The bounds of the control input are very relevant!! If too big, the COM
            %allows it but the SOM not (the COM is a spring that from a time step
            %to the next one can be mathematically infinetely stretched). The
            %bounds should be approximately the maximum displacement expected.
            lbw = [lbw; -ubound*ones(n_u * M,1)];
            ubw = [ubw;  ubound*ones(n_u * M,1)];

            delta_u = [U(:, 1) - U0, diff(U, 1, 2)];
            delta_u = delta_u;
            delta_u = delta_u(:);

            node_l_delta = U([1, 3, 5], :);
            node_r_delta = U([2, 4, 6], :);
            node_l = [node_l_delta(:, 1) + U0([1, 3, 5])];
            node_r = [node_r_delta(:, 1) + U0([2, 4, 6])];
            for k = 2:M
                node_l = [node_l, node_l_delta(:, k) + node_l(:, k - 1)];
                node_r = [node_r, node_r_delta(:, k) + node_r(:, k - 1)];
            end
            for k = 1:M
                g = [g; sum((node_l(:, k) - node_r(:, k)).^2) - lCloth^2];
                gbound = 0;
                lbg = [lbg; -gbound];
                ubg = [ubg;  gbound];
            end
            U = U';
            U = U(:);
            w = [w; U];
            prediction = P(:, 1) + S_tilde * U + D_tilde * P(1:n_d * H, 3);
            for k = (1:H)
                g = [g;
                     Cc([k, H + k, 2*H + k, 3*H + k, 4*H + k, 5*H + k], [k, H + k, 2*H + k, 3*H + k, 4*H + k, 5*H + k]) *...
                     prediction([k, H + k, 2*H + k, 3*H + k, 4*H + k, 5*H + k]) - bc([k, H + k, 2*H + k, 3*H + k, 4*H + k, 5*H + k])];
                ubg = [ubg;
                       zeros(n_y, 1)];
                lbg = [lbg;
                       -inf * ones(n_y, 1)];
            end
            error = prediction - P(:, 2); % F_tilde + S_tilde * delta_u - reference + D_tilde * delta
            % Update objective function
            obj = obj + error'*Q*error;
            obj = obj + R*(U'*U);

            qp_prob = struct('f', obj, 'x', w, 'g', g, 'p', P);
            opts = struct;
            opts.print_time = 0;
            opts.ipopt.print_level = 0; %0 to print the minimum, 3 to print the maximum
            opts.ipopt.warm_start_init_point = 'yes'; %warm start

            solver = nlpsol('solver', 'ipopt', qp_prob, opts);

            %% Run MPC simulation loop


            % plot(phi_l_Traj);
            u_ini = x_ini_SOM(SOM.coord_ctrl);
            u_bef = u_ini;

            % Initialize things
            reference = zeros(H * n_y, 1);
            store_state(:,end + 1) = x_ini_SOM;


            tT = 0;

            time_vector = (1:size(phi_l_Traj,1) + 1)' * COM.dt;
            e = exp(1);
            wind_velocity = 2.4 * (e * ones(length(time_vector), 1)).^(-((time_vector - 5) ./ 2).^2);
            if disturbance
                wind_force = -0.5 * 1.225 * wind_velocity.^2 * 0.09; 
                noise_var = 0.005;
                sigma_epsilon_i = noise_var*eye(length(wind_force));
                padded_identity = eye(H, length(wind_force));
                Nu = circshift(padded_identity, [0, 1]);
                Xi = circshift(padded_identity, [0, 0]);
                sigma_nu_i = (Nu - Xi) * sigma_epsilon_i * (Nu - Xi)';
                sigma_nu_tilde = blkdiag(sigma_nu_i, sigma_nu_i, sigma_nu_i);
                sigma_y_tilde = D_tilde * sigma_nu_tilde * D_tilde';
            else
                wind_force = 0 * wind_velocity; % acts along the y axis
                sigma_tilde_hy = zeros(n_y * H);
            end
            writematrix([zeros(length(wind_force), 1), wind_force, zeros(length(wind_force), 1)], "disturbance_profile.csv");
            % plot(time_vector, wind_force)
            xnew = x_ini_COM_rot;
            free_response = [x_ini_COM_rot'];
            for k = (2: 10*N)
                xnew = A * xnew + COM.dt .* ext_force;
                %free_response(k, :) = xnew';
            end
            for k = (2: N)
                xnew = A * xnew + COM.dt .* ext_force;
                free_response(k, :) = xnew';
            end
    %         plot(free_response);
    %         ficurr = x_ini_COM_rot(coord_l)';
    %         ficurr = repmat(ficurr, [N, 1]);
            ficurr = free_response(1:N, coord_l);
            ficurr = ficurr(:);
            last_disturbance = zeros(n_d * H, 1);
            last_disturbance(H + 1) = wind_force(1);
            last_disturbance(H + 2) = wind_force(2) - wind_force(1);
            ficurr = ficurr + D_tilde_steady_state * last_disturbance;
            d = zeros(H * n_y, 1);
            f_tilde_i = Psi * ficurr + d;
            SOM_gravity = SOM.Fg;
            delta_u_bef = zeros(6, 1);
            u_bef_new_frame = Rcloth^-1 * [u_bef([1, 3, 5]), u_bef([2, 4, 6])];
            u_bef_new_frame = reshape(u_bef_new_frame', 1, 6)';
            for i=2:size(phi_l_Traj,1)
                if i>=size(phi_l_Traj,1)-H-1 %the last N timesteps, trajectory should remain constant
                    Traj_l_Hp = repmat(phi_l_Traj(end,:),H, 1);
                    Traj_r_Hp = repmat(phi_r_Traj(end,:),H, 1);
                else
                    Traj_l_Hp = phi_l_Traj(i:i+H-1,:);
                    Traj_r_Hp = phi_r_Traj(i:i+H-1,:);
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
    %             plot(Traj_l_Hp_rot)
    %             hold on
    %                         plot(Traj_r_Hp_rot)
    % hold off
                reference(0 * H + 1: 1 * H, :) = Traj_l_Hp_rot(:,1)';
                reference(1 * H + 1: 2 * H, :) = Traj_r_Hp_rot(:,1)';
                reference(2 * H + 1: 3 * H, :) = Traj_l_Hp_rot(:,2)';
                reference(3 * H + 1: 4 * H, :) = Traj_r_Hp_rot(:,2)';
                reference(4 * H + 1: 5 * H, :) = Traj_l_Hp_rot(:,3)';
                reference(5 * H + 1: 6 * H, :) = Traj_r_Hp_rot(:,3)';
                sliced_wind_force = wind_force(i-1:min(i+H-1, end));

                if i == 3
                    sliced_wind_force(1) = 0;
                end
                delta_sliced_wind_force = diff(sliced_wind_force);
                wind_force_x = -1/3 * delta_sliced_wind_force;
                wind_force_z = zeros(length(delta_sliced_wind_force), 1);


                wind_force_x = [wind_force_x;
                                zeros(max(0, H - length(wind_force_x)), 1)];
                delta_sliced_wind_force = [delta_sliced_wind_force;
                                           zeros(max(0, H - length(delta_sliced_wind_force)), 1)];
                wind_force_z = [wind_force_z;
                                zeros(max(0, H - length(wind_force_z)), 1)];

                predicted_disturbance = [wind_force_x;
                                         delta_sliced_wind_force;
                                         wind_force_z];
                predicted_disturbance = [predicted_disturbance;
                                         zeros((n_y - n_d) * H, 1)];

                padded_ubef = [u_bef_new_frame; zeros(length(f_tilde_i) - length(delta_u_bef), 1)];

                % Initial seed of the optimization (for "u")
                args.x0  = zeros(6*M, 1);
                % Build constraint matrix and vector.

                if i < t1 - 1
                    Cconstr = [[eye(2*H, 2*H), zeros(2*H, 4*H)];
                               zeros(4*H, 6*H)];
                    bconstr = [repmat(0.35, 2*H, 1);
                               zeros(4*H, 1)];
                elseif i >= t1 - 1 && i <= t1 + t2 + t3 + 1
                    Cconstr = [zeros(2*H, 6*H);
                               [zeros(2*H, 2*H), -1*eye(2*H, 2*H), zeros(2*H, 2*H)];
                               zeros(2*H, 6*H)];
                    bconstr = [zeros(2*H, 1);
                               repmat(-0.25, 2*H, 1);
                               zeros(2*H, 1)];
                elseif i > t1 + t2 + t3 + 1
                    Cconstr = [[-1*eye(2*H, 2*H), zeros(2*H, 4*H)];
                               zeros(4*H, 6*H)];
                    bconstr = [repmat(-0.55, 2*H, 1);
                               zeros(4*H, 1)];
                end
                sigma_tilde_hy = Cconstr * sigma_y_tilde * Cconstr';
                if chance
                    correction = norminv(1 - alpha_err/(2 * n_y * H)) * diag(sigma_tilde_hy).^0.5;
                    bconstr = bconstr - correction;
                end
                args.p = [f_tilde_i, reference, predicted_disturbance, padded_ubef, sparse(Cconstr), sparse(bconstr)];
                t0 = tic;
                sol = solver('x0', args.x0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg, 'p', args.p);
                tmp = full(sol.x);
                stats = solver.stats();
                if stats.success == 0
                    "Optimization failed"
                    return
                end
                sol = full(sol.x);
                u_rot =  reshape(sol,M,6)';
                u_sliced = u_rot(:, 1);            

                delta_u_bef = u_sliced; % Optimization variable is displacement delta_u; then we penalize the upper corner accel.
                u_rot_l = u_sliced([1 3 5], 1);
                u_rot_r = u_sliced([2 4 6], 1);
                u_base = Rcloth * [u_rot_l, u_rot_r];
                u = reshape(u_base',[6,1]);
                u_SOM = u_bef + u;% + [0.0005; 0.0005; -0.0005; -0.0005; 0.0001; 0.0001];

                d_debug = u_SOM - u_bef;

                delta_u_chosen_free = zeros(n_u * M, 1);
    %             d_debug = Rcloth^-1 * [d_debug([1, 3, 5]), d_debug([2, 4, 6])];
    %             d_debug = reshape(d_debug', 1, 6)';
                for j = (1:n_u)
                    delta_u_chosen_free((j - 1) * M + 1) = u_sliced(j);
                end
                last_disturbance = zeros(n_d * H, 1);
                last_disturbance(1) = predicted_disturbance(1,  1);
                last_disturbance(H + 1) = predicted_disturbance(2, 1);

                % Update free step response by applying the first input only
                ficurr = Gamma * ficurr + S_tilde_steady_state * delta_u_chosen_free...
                        + D_tilde_steady_state * last_disturbance;  
                ficurr_rsh = reshape(ficurr', N, 6)';
                ficurr_base = Rcloth * [ficurr_rsh([1, 3, 5], :), ficurr_rsh([2, 4, 6], :)];
                ficurr_base = [ficurr_base(1, 1:N)';
                               ficurr_base(1, N+1:end)';
                               ficurr_base(2, 1:N)';
                               ficurr_base(2, N+1:end)';
                               ficurr_base(3, 1:N)';
                               ficurr_base(3, N+1:end)';];

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
                next_state_SOM = [next_pos_SOM; next_vel_SOM];
                phi_ini_SOM = full(next_pos_SOM);
                dphi_ini_SOM = full(next_vel_SOM);  

                measured_output = next_state_SOM(coord_nl);
                cloth_x = u_SOM([2 4 6]) - u_SOM([1 3 5]);
                cloth_y = [-cloth_x(2) cloth_x(1) 0]';
                cloth_z = cross(cloth_x,cloth_y);

                cloth_x = cloth_x/norm(cloth_x);
                cloth_y = cloth_y/norm(cloth_y);
                cloth_z = cloth_z/norm(cloth_z);
                Rcloth = [cloth_x cloth_y cloth_z];
                u_bef_new_frame = Rcloth^-1 * [u_bef([1, 3, 5]), u_bef([2, 4, 6])];
                u_bef_new_frame = reshape(u_bef_new_frame', 1, 6)';
                phi_ini_SOM = full(next_state_SOM(1:3*SOM.n_nodos));
                dphi_ini_SOM = full(next_state_SOM((1+3*SOM.n_nodos):6*SOM.n_nodos));  

                % Compute the model mismatch.
                d = [];
                pred_output = [];
                for j = (1:n_y)
                    fkk = ficurr_base((j-1) * N + 1);
                    pred_output = [pred_output, ficurr_base((j-1) * N + 1:j * N)];
                    curr_output_channel = measured_output(j);
                    d = [d; repmat(curr_output_channel - fkk, [H, 1])];
                end
                left = pred_output(:, [1, 3, 5]);
                right = pred_output(:, [2, 4, 6]);
    %             scatter3(Traj_l_Hp(:,1)', Traj_l_Hp(:,2)', Traj_l_Hp(:,3)', 'LineWidth',1.5);
    %             hold on
    %             scatter3(Traj_r_Hp(:,1)', Traj_r_Hp(:,2)', Traj_r_Hp(:,3)', 'LineWidth',1.5);
%                 plot3(left(:, 1), left(:, 2), left(:, 3));
%                 hold on
%                 plot3(right(:, 1), right(:, 2), right(:, 3));
%                 scatter3(measured_output(1), measured_output(3), measured_output(5), 'LineWidth',1.5);
%                 scatter3(measured_output(2), measured_output(4), measured_output(6), 'LineWidth',1.5);

%                 scatter3(u_SOM_benchmark(1), u_SOM_benchmark(3), u_SOM_benchmark(5), 'LineWidth',1.5);
%                 scatter3(u_SOM_benchmark(2), u_SOM_benchmark(4), u_SOM_benchmark(6), 'LineWidth',1.5);
%     
%                 scatter3(u_SOM(1), u_SOM(3), u_SOM(5), 'LineWidth',1.5);
%                 scatter3(u_SOM(2), u_SOM(4), u_SOM(6), 'LineWidth',1.5);
                f_tilde_i_base = Psi * ficurr_base + d;
                ftilde_rsh = reshape(f_tilde_i_base', H, 6)';
                f_tilde_i = (Rcloth^-1) * [ftilde_rsh([1, 3, 5], :), ftilde_rsh([2, 4, 6], :)];
                f_tilde_i = [f_tilde_i(1, 1:H)';
                               f_tilde_i(1, H+1:end)';
                               f_tilde_i(2, 1:H)';
                               f_tilde_i(2, H+1:end)';
                               f_tilde_i(3, 1:H)';
                               f_tilde_i(3, H+1:end)';];
                ficurr_base_rsh = reshape(ficurr_base', N, 6)';
                ficurr = (Rcloth^-1) * [ficurr_base_rsh([1, 3, 5], :), ficurr_base_rsh([2, 4, 6], :)];
                ficurr = [ficurr(1, 1:N)';
                           ficurr(1, N+1:end)';
                           ficurr(2, 1:N)';
                           ficurr(2, N+1:end)';
                           ficurr(3, 1:N)';
                           ficurr(3, N+1:end)';];
                [reduced_pos,~] = take_reduced_mesh(phi_ini_SOM, ...
                                    dphi_ini_SOM, ...
                                    nx, COM.row);
%                 scatter3(reduced_pos(1:nr*nc), reduced_pos(nr*nc + 1:nr*nc*2), reduced_pos(nr*nc*2 + 1:nr*nc*3));
%                 ylim([-0.4, 0.7]);
%                 xlim([-0.4, 1.2]);
%                 zlim([-0.2, 0.3]);
%                 hold off
                % Store things
                store_state(:,end + 1) = [phi_ini_SOM;dphi_ini_SOM];
    %             if i > 100
    %                 store_state = lowpass(store_state',1, 50, ImpulseResponse='iir')';
    %                 %figure
    %                 lowpass(store_state(1:3, :)',1, 50, ImpulseResponse='iir')
    %             end
                store_u(:,end + 1) = u_SOM;    

                if(mod(i,100)==0)
                    fprintf(['Iter: ', num2str(i), ...
                             ' \t Avg. time/iter: ', num2str(tT/i*1000), ' ms \n']);
                end
            end

            tT = tT+toc(t0);
            time = 0:0.01:size(store_state,2)*0.01-0.01;
            plot(store_state(coord_nl(2),:)', store_state(coord_nl(4),:)','linewidth',2.5);
            hold on
            ylim([-0.05, 0.4]);

            %plot(store_state(coord_nl(2),:)',store_state(coord_nl(4),:)','linewidth',2.5);


    %         fig1 = figure(1);
    %         fig1.Color = [1,1,1];
    %         fig1.Units = 'normalized';
    %         fig1.Position = [0.5 0.6 0.5 0.3];
    % 
    %         subplot(7,2,1:2:12);
    %         pa1som=plot(time,store_state(coord_nl([1 3 5]),:)','linewidth',2.5);
    %         hold on
    %         pa1ref=plot(time,phi_l_Traj(:,:),'--k','linewidth',2.2);
    %         hold off
    %         title('\textbf{Lower Left Corner}', 'Interpreter', 'latex')
    %         grid on
    %         xlabel('Time [s]', 'Interpreter', 'latex')
    %         ylabel('Position [m]', 'Interpreter', 'latex')
    %         xlim([0 round(time(end))])
    %         set(gca, 'TickLabelInterpreter', 'latex');
    % 
    %         subplot(7,2,2:2:12);
    %         plot(time,store_state(coord_nl([2 4 6]),:)','linewidth',2.5)
    %         hold on
    %         plot(time,phi_r_Traj(:,:)','--k','linewidth',2.2)
    %         hold off
    %         title('\textbf{Lower Right Corner}', 'Interpreter', 'latex')
    %         grid on
    %         xlabel('Time [s]', 'Interpreter', 'latex')
    %         ylabel('Position [m]', 'Interpreter', 'latex')
    %         xlim([0 round(time(end))])
    %         set(gca, 'TickLabelInterpreter', 'latex');
    % 
    %         Lgnd1 = legend([pa1som' pa1ref(1)], ...
    %                        '$p_x\ (y_1,y_2)$','$p_y\ (y_3,y_4)$', '$p_z\ (y_5,y_6)$', '$r$', ...
    %                        'Orientation','horizontal', 'Interpreter', 'latex');
    %         Lgnd1.Position(1) = 0.5-Lgnd1.Position(3)/2;
    %         Lgnd1.Position(2) = 0.02;
    %         Lgnd1.Position(3:4) = Lgnd1.Position(3:4) + 0.01;

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
            runtimes(end+1) = tT/size(phi_l_Traj,1)*1000;
            errors(end+1) = 1000*eMAEm;
            tottimes(end+1) = tT;
        end
        tot_runtimes = [tot_runtimes; runtimes];
        tot_errors = [tot_errors; errors];
        if chance == false
            break
        end
    end
end
plot(phi_r_Traj(:,1)',phi_r_Traj(:,2)','--k','linewidth',2.2)
%plot(phi_l_Traj(:,1)',phi_l_Traj(:,2)','--k','linewidth',2.2)
line([0.35, 0.55], [0.25, 0.25], 'LineWidth', 1.5, 'Color', 'blue');
line([0.55, 0.55], [0.25, 0.05], 'LineWidth', 1.5, 'Color', 'blue');
line([0.55, 0.35], [0.05, 0.05], 'LineWidth', 1.5, 'Color', 'blue');
line([0.35, 0.35], [0.05, 0.25], 'LineWidth', 1.5, 'Color', 'blue');

x = [0.35 0.55 0.55 0.35];
y = [0.25 0.25 0.05 0.05];
fill(x,y,'b', 'FaceAlpha', alpha)

legend('$\alpha = 1$', '$\alpha = 0.5$', '$\alpha = 0.1$', 'No ch.\ constr.', '$r$', 'Obst.', 'Interpreter', 'latex', 'Location', 'north outside', 'NumColumns', 6);
set(gca, 'TickLabelInterpreter', 'latex');
xlabel('$x$ [m]', 'interpreter', 'latex');
ylabel('$y$ [m]', 'interpreter', 'latex');
ylim([-0.05, 0.4]);
grid on
hold off

% writematrix(tot_errors, "obstacle_qdmc_error_sim_traj" + ".csv");
% writematrix(tot_runtimes, "obstacle_qdmc_time_sim_traj" +".csv");

    
% writematrix(error_trajs, "all_trajs_error.csv");
% writematrix(runtime_trajs, "all_trajs_opti_time.csv")