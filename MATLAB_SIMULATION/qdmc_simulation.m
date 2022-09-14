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
seed = 1; % Set it to 1 in the deterministic version
error_seeds = [];

rng(seed);
disturbance = false;
% Define COM parameters (computation oriented model)
COM = struct;
COM.row = 4; COM.col = 4;
COM.mass = 0.1;
COM.grav = 9.8;
COM.dt = 0.02;

% parameters from calibration    
COM.stiffness = [-41.7998405537012  -39.1879707090799 -91.4967890838762]; 
COM.damping = [-1.72864429464933   -1.33085664353071   -2.2554123678956]; 
COM.z_sum = 0.0709542759848484;

% Define the SOM
nx = 10; ny = nx;
[SOM,pos] = initialize_model_parameters(nx,ny);

% Definite initial position of the nodes (defined here because we need it
% to compute ext_force
x_ini_SOM = [reshape(pos,[3*nx*ny 1]);zeros(3*nx*ny,1)]; %initial velocity=0
[reduced_pos,~] = take_reduced_mesh(x_ini_SOM(1:300),x_ini_SOM(301:600));
x_ini_COM = [reduced_pos;zeros(48,1)];

% Initial position of the COM nodes
COM.nodeInitial = lift_z(reshape(x_ini_COM(1:48), [16,3]), COM);

% Find initial spring length in each direction x,y,z
[COM.mat_x, COM.mat_y, COM.mat_z] = compute_l0_linear(COM,0);

% Find linear matrices
[A, B, ext_force] = create_model_linear_matrices(COM);

% Build disturbance matrix (wind is 3-dimensional disturbance)
B_delta = zeros(2 * COM.row * COM.col * 3, 3);
B_delta(COM.row * COM.col * 3 + 1:COM.row * COM.col * 4, 1) = COM.dt * 1/COM.mass * ones(COM.row * COM.col, 1);
B_delta(COM.row * COM.col * 4 + 1:COM.row * COM.col * 5, 2) = COM.dt * 1/COM.mass * ones(COM.row * COM.col, 1);
B_delta(COM.row * COM.col * 5 + 1:COM.row * COM.col * 6, 3) = COM.dt * 1/COM.mass * ones(COM.row * COM.col, 1);
B_delta([13+48 16+48], 1) = zeros(2,1);
B_delta([29+48 32+48], 2) = zeros(2,1);
B_delta([45+48 48+48], 3) = zeros(2,1);
B_delta = sparse(B_delta);
nr = COM.row;
nc = COM.col;

coord_l = [1 nc 1+nr*nc nr*nc+nc 2*nr*nc+1 2*nr*nc+nc]; %coordinates of the linear model (of the lower corners)
coord_nl = [1 nx 1+nx*ny nx*ny+nx 2*nx*ny+1 2*nx*ny+nx]; %coordinates of the non-linear model (of the lower corners)

%% Step response for FSR model

% Compute the coefficients of the step response of MIMO system.
% Note that the system considered is linear, thus it suffices to 
% apply one step for each input-output pair.
N = 300; % finite-step-response horizon until steady state is reached.
H = 30;
Ms = (2:30); % Control horizon, set it to 5 to avoid the tuning of this parameter
runtimes = [];
errors = [];
tottimes = [];
for M = Ms
    n_u = 6;
    n_y = 6;
    step_amplitude = 1e10; % Offset w.r.t. initial position
    step_response = zeros(N, length(x_ini_COM));
    free_response = zeros(10 * N, length(x_ini_COM));

    MIMO_step_response = zeros(0, N, n_y); % MIMO_step_response(a, b, c) is the output response of channel c at timestep b
                                           % given a step input in channel a

    % WE MUST START AT THE EQUILIBRIUM
    xnew = x_ini_COM;
    %xnew = A * xnew + B * [0.05; 0; 0; -0.5; 0; 0] + COM.dt .* ext_force;
    for k = (1: 10 * N)
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
    end

    % Now, identify the transfer from each disturbance channel (3) to each
    % output channel
    %figure
    n_d = 3;
    step_amplitude = 1e-5; % Offset w.r.t. initial position
    step_response = zeros(N, length(x_ini_COM));
    free_response = zeros(10 * N, length(x_ini_COM));

    MIMO_step_response_disturbance = zeros(0, N, n_y); % MIMO_step_response(a, b, c) is the output response of channel c at timestep b
                                           % given a step input in channel a

    % WE MUST START AT THE EQUILIBRIUM
    xnew = x_ini_COM;
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
    S_tilde(S_tilde<1e-12) = 0;
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
    P = SX.sym('fk', H * n_y, 3); % Parameters. First col is f_tilde_k, second 
                                  % col is reference signal, third col is
                                  % predicted disturbance
    ubound = 5e-3; %0.8e-3, 50e-3
    % Tune R and Q; normalize w.r.t. maximum input variation and maximum est.
    % error
    R = 1e7; % Control variation in the order of [mm] (favor small variations - 
             % the simulator is a bottleneck! Also, too harsh movements of the cloth
             % seem to cause oscillations).
    Q = 1e4; % Tracking error in the order of [mm] (from LMPC).

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
    T = 699;
    % New NLP variable for the control
    Uk = SX.sym('U', n_u * M);
    w = [w; Uk];
    %The bounds of the control input are very relevant!! If too big, the COM
    %allows it but the SOM not (the COM is a spring that from a time step
    %to the next one can be mathematically infinetely stretched). The
    %bounds should be approximately the maximum displacement expected.
    lbw = [lbw; -ubound*ones(n_u * M,1)];
    ubw = [ubw;  ubound*ones(n_u * M,1)];

    error = P(:, 1) + S_tilde * Uk - P(:, 2) + D_tilde * P(1:n_d * H, 3); % F_tilde + S_tilde * delta_u - reference + D_tilde * delta
    % Update objective function

    obj = obj + error'*Q*error;
    obj = obj + R*(Uk'*Uk);

    qp_prob = struct('f', obj, 'x', w, 'g', g, 'p', P);
    opts = struct;
    opts.print_time = 0;
    opts.ipopt.print_level = 0; %0 to print the minimum, 3 to print the maximum
    opts.ipopt.warm_start_init_point = 'yes'; %warm start

    solver = nlpsol('solver', 'ipopt', qp_prob, opts);

    %% Run MPC simulation loop

    phi_l_Traj = load('required_files/trajectories/phi_l_3D.csv');
    phi_r_Traj = load('required_files/trajectories/phi_r_3D.csv');
    lifting = false; %Set to true to test the trajecotry in our simulation experiments; false to tune the horizons
    if lifting
        phi_l_Traj(:, 2) = phi_l_Traj(1, 2) + (1:length(phi_l_Traj)).^2 / length(phi_l_Traj)^2 * 0.2;
        phi_r_Traj(:, 2) = phi_r_Traj(1, 2) + (1:length(phi_r_Traj)).^2/ length(phi_r_Traj)^2 * 0.2;
        phi_l_Traj(:, 3) = phi_l_Traj(1, 3) + (1:length(phi_l_Traj)).^2 / length(phi_l_Traj)^2 * 0.2;
        phi_r_Traj(:, 3) = phi_r_Traj(1, 3) + (1:length(phi_r_Traj)).^2/ length(phi_r_Traj)^2 * 0.2;

        phi_l_Traj(500:end, 2) = phi_l_Traj(500, 2);
        phi_r_Traj(500:end, 2) = phi_l_Traj(500, 2);
        phi_l_Traj(500:end, 3) = phi_l_Traj(500, 3);
        phi_r_Traj(500:end, 3) = phi_l_Traj(500, 3);
    end
    % plot(phi_l_Traj);
    u_ini = x_ini_SOM(SOM.coord_controlados);
    u_bef = u_ini;


    % Counters for the loop
    sHp = 1; % start prediction horizon
    Hp = H; % end prediction horizon

    % Initialize things
    reference = zeros(H * n_y, 1);
    store_state(:,1) = x_ini_SOM;
    store_state(:,2) = x_ini_SOM;


    tT = 0;

    time_vector = (1:size(phi_l_Traj,1) + 1)' * COM.dt;
    e = exp(1);
    wind_velocity = 5 / 3.6 * (e * ones(length(time_vector), 1)).^(-((time_vector - 3) ./ 1^2).^2);
    if disturbance
        wind_force = 0.5 * 1.225 * wind_velocity.^2 * 0.09; % acts along the y axis
    else
        wind_force = 0 * wind_velocity; % acts along the y axis
    end

    writematrix([zeros(length(wind_force), 1), wind_force, zeros(length(wind_force), 1)], "disturbance_profile.csv");
%     plot(time_vector, wind_force)
    ficurr = x_ini_COM(coord_l)';
    ficurr = repmat(ficurr, [N, 1]);
    last_disturbance = zeros(n_d * H, 1);
    last_disturbance(H + 1) = wind_force(1);
    last_disturbance(H + 2) = wind_force(2) - wind_force(1);
    ficurr = ficurr(:) + D_tilde_steady_state * last_disturbance;
    d = zeros(H * n_y, 1);
    f_tilde_i = Psi * ficurr + d;
    measured_output = x_ini_SOM(coord_nl);
    SOM_gravity = SOM.Fg;
    for i=3:size(phi_l_Traj,1)
        Hp = Hp+1;
        sHp = sHp+1;
        if i>=size(phi_l_Traj,1)-H % in the last H timesteps, trajectory should 
                                   % remain constant
            Traj_l_Hp = repmat(phi_l_Traj(end,:),H);
            Traj_r_Hp = repmat(phi_r_Traj(end,:),H);
        else
            Traj_l_Hp = phi_l_Traj(sHp:Hp,:);
            Traj_r_Hp = phi_r_Traj(sHp:Hp,:);
        end

        % Define reference in the prediction horizon (moving window)
        reference(0 * H + 1: 1 * H, :) = Traj_l_Hp(:,1)';
        reference(1 * H + 1: 2 * H, :) = Traj_r_Hp(:,1)';
        reference(2 * H + 1: 3 * H, :) = Traj_l_Hp(:,2)';
        reference(3 * H + 1: 4 * H, :) = Traj_r_Hp(:,2)';
        reference(4 * H + 1: 5 * H, :) = Traj_l_Hp(:,3)';
        reference(5 * H + 1: 6 * H, :) = Traj_r_Hp(:,3)';
        sliced_wind_force = wind_force(i-1:min(i+H-1, end));
%         plot(sliced_wind_force);
%         hold on
        if i == 3
            sliced_wind_force(1) = 0;
        end
        delta_sliced_wind_force = diff(sliced_wind_force);
        wind_force_x = zeros(length(delta_sliced_wind_force), 1);
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


        % Initial seed of the optimization (for "u")
        args.x0  = [repmat([zeros(n_u,1)], M,1)];
        args.p = [f_tilde_i, reference, predicted_disturbance];
        t0 = tic;
        sol = solver('x0', args.x0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg, 'p', args.p);
        stats = solver.stats();
        if stats.success == 0
            "Optimization failed"
            return
        end
        u = full(sol.x(:, :)); 
        delta_u_chosen = u([1, M + 1, 2*M + 1, 3*M + 1, 4*M + 1, 5*M + 1]);
        u_SOM = delta_u_chosen + u_bef;
        u_bef = u_SOM; %update
        tT=tT+toc(t0);
        SOM_perturbation = [zeros(nx * ny, 1);
                           repmat(sliced_wind_force(1), nx*ny, 1);
                           zeros(nx * ny, 1)];
        %hold on
        SOM_perturbation = SOM_perturbation / COM.mass;
        SOM_perturbation = SOM.Mlum*SOM_perturbation*SOM.rho;
        SOM.Fg = SOM_gravity + SOM_perturbation; % Wind acting on the y axis.

        [next_state_SOM] = cloth_simulator_secondorder([store_state(:,i-1);store_state(:,i-2)],u_SOM,SOM);
        measured_output = next_state_SOM(coord_nl);
        phi_ini_SOM = full(next_state_SOM(1:3*SOM.n_nodos));
        dphi_ini_SOM = full(next_state_SOM((1+3*SOM.n_nodos):6*SOM.n_nodos));    

        delta_u_chosen_free = zeros(n_u * M, 1);
        for j = (1:n_u)
            delta_u_chosen_free((j - 1) * M + 1) = delta_u_chosen(j);
        end
        last_disturbance = zeros(n_d * H, 1);
        last_disturbance(H + 1) = delta_sliced_wind_force(1);

        % Update free step response by applying the first input only
        ficurr = Gamma * ficurr + S_tilde_steady_state * delta_u_chosen_free...
                + D_tilde_steady_state * last_disturbance;  
        % Compute the model mismatch.
        d = [];
        for j = (1:n_y)
            fkk = ficurr((j-1) * N + 1);
            curr_output_channel = measured_output(j);
            d = [d; repmat(curr_output_channel - fkk, [H, 1])];
        end
        f_tilde_i = Psi * ficurr + d;


        % Store things
        store_state(:,i) = [phi_ini_SOM;dphi_ini_SOM];
        store_u(:,i) = u_SOM;    

        if(mod(i,100)==0)
            fprintf(['Iter: ', num2str(i), ...
                     ' \t Avg. time/iter: ', num2str(tT/i*1000), ' ms \n']);
        end
    end

    tT = tT+toc(t0);
%         time = 0:0.01:size(store_state,2)*0.01-0.01;
% 
%     fig1 = figure(1);
%     fig1.Color = [1,1,1];
%     fig1.Units = 'normalized';
%     fig1.Position = [0.5 0.6 0.5 0.3];
% 
%     subplot(7,2,1:2:12);
%     pa1som=plot(time,store_state(coord_nl([1 3 5]),:)','linewidth',1.5);
%     hold on
%     pa1ref=plot(time,phi_l_Traj(:,:),'--k','linewidth',1.2);
%     hold off
%     title('\textbf{Lower Left Corner}', 'Interpreter', 'latex')
%     grid on
%     xlabel('Time [s]', 'Interpreter', 'latex')
%     ylabel('Position [m]', 'Interpreter', 'latex')
%     xlim([0 round(time(end))])
%     set(gca, 'TickLabelInterpreter', 'latex');
% 
%     subplot(7,2,2:2:12);
%     plot(time,store_state(coord_nl([2 4 6]),:)','linewidth',1.5)
%     hold on
%     plot(time,phi_r_Traj(:,:)','--k','linewidth',1.2)
%     hold off
%     title('\textbf{Lower Right Corner}', 'Interpreter', 'latex')
%     grid on
%     xlabel('Time [s]', 'Interpreter', 'latex')
%     ylabel('Position [m]', 'Interpreter', 'latex')
%     xlim([0 round(time(end))])
%     set(gca, 'TickLabelInterpreter', 'latex');
% 
%     Lgnd1 = legend([pa1som' pa1ref(1)], ...
%                    '$p_x\ (y_1,y_2)$','$p_y\ (y_3,y_4)$', '$p_z\ (y_5,y_6)$', '$r$', ...
%                    'Orientation','horizontal', 'Interpreter', 'latex');
%     Lgnd1.Position(1) = 0.5-Lgnd1.Position(3)/2;
%     Lgnd1.Position(2) = 0.02;
%     Lgnd1.Position(3:4) = Lgnd1.Position(3:4) + 0.01;

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
% writematrix(errors, "qdmc_error_sim.csv");
% writematrix(runtimes, "qdmc_time_sim.csv");
% writematrix(tottimes, "qdmc_tot_time_sim.csv");