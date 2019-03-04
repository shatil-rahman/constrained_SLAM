%Batch LQG Estimation Exercise 2
%Author: Shatil Rahman

clear all
close all
clc

%% Initializations
no_experiments = 50;
MC_errors = zeros(no_experiments,4);

for counter=1:no_experiments
    load dataset2.mat
    disp(['Experiment no ',num2str(counter)])
    
    row = 1;
    k_start = 4000;
    k_step = 2;
    k_end = 5000;

    X_0 = [x_true(k_start); y_true(k_start); wrapToPi(th_true(k_start))]; 
    P_0 = [0.01, 0.0, 0.0;
           0.0, 0.01, 0.0
           0.0, 0.0, 1.0];  

    X_L_0 = 3*randn(2*length(l),1);
    for i=1:length(l)
        X_L_0(row) = l(i,1) + 0.25*randn();
        X_L_0(row+1) = l(i,2) + 0.25*randn();
        row = row + 2;
    end

    %% Noise levels
    tic
    % Odometry Noise
    noise_odom = 0.0;%0.02;
    noise_v = 0.0;%0.05;

    om(k_start:k_step:k_end) = om(k_start:k_step:k_end) + noise_odom*randn(length(om(k_start:k_step:k_end)),1);
    v(k_start:k_step:k_end) = v(k_start:k_step:k_end) + noise_v*randn(length(v(k_start:k_step:k_end)),1);

    % Range/bearing noise
    noise_r = 0.0;%0.8;
    noise_b = 0.0;%0.1;

    for i=k_start:k_step:k_end
        for j = 1:17
            if r(i,j) ~=0
                r(i,j) = r(i,j) + noise_r*randn();
                b(i,j) = b(i,j) + noise_b*randn();
            end
        end
    end

    Q_estimator = [v_var + noise_v, 0;
                   0, om_var + noise_odom];
    R_estimator = [r_var + noise_r, 0;
                    0, b_var + noise_b];
    % constrained_landmarks = [17,2;  2,16; 7,10; 17,16; 10,3; 1,2; 3,17];

    %Constraint info

    constr_noise = 0.0;
    [constrained_landmarks, fixed_landmarks] = edge_counter(l,r(k_start:k_step:k_end,:));
    Y_const = zeros(length(constrained_landmarks),1);
    for i = 1:length(constrained_landmarks)
       Y_const(i) = norm(l(constrained_landmarks(i,1),:) -  l(constrained_landmarks(i,2),:)) + constr_noise*randn();
    end

    %%
    %tic
    [q_batch,L_batch, del_Xs_batch] = GN_Estimator( X_0, P_0,X_L_0, Q_estimator, R_estimator, v(k_start:k_step:k_end),...
         om(k_start:k_step:k_end), r(k_start:k_step:k_end,:),b(k_start:k_step:k_end,:), t(k_start:k_step:k_end),d,...
         x_true(k_start:k_step:k_end), y_true(k_start:k_step:k_end), th_true(k_start:k_step:k_end),l);
    %toc
    %%
    %tic
    [q_prob,L_prob, del_Xs_prob] = prob_SLAM( X_0, P_0,X_L_0, Q_estimator, R_estimator, v(k_start:k_step:k_end),...
        om(k_start:k_step:k_end), r(k_start:k_step:k_end,:),b(k_start:k_step:k_end,:),Y_const, t(k_start:k_step:k_end),d,...
           constrained_landmarks,...
         x_true(k_start:k_step:k_end), y_true(k_start:k_step:k_end), th_true(k_start:k_step:k_end),l);
    %toc
    %% 
    %tic
    [q_eq,L_eq, del_Xs_eq] = constrained_SLAM( X_0, P_0,X_L_0, Q_estimator, R_estimator, v(k_start:k_step:k_end),...
        om(k_start:k_step:k_end), r(k_start:k_step:k_end,:),b(k_start:k_step:k_end,:),Y_const, t(k_start:k_step:k_end),d,...
           constrained_landmarks,...
         x_true(k_start:k_step:k_end), y_true(k_start:k_step:k_end), th_true(k_start:k_step:k_end),l);
    %toc
    %%
    %tic
    X_L_0_fixed = X_L_0;
    % X_L_0_fixed(2*fixed_landmarks-1) = l(fixed_landmarks,1);
    % X_L_0_fixed(2*fixed_landmarks) = l(fixed_landmarks,2);

    [q_var_rem,L_var_rem, del_Xs_var_rem] = var_rem_SLAM( X_0, P_0,X_L_0_fixed, Q_estimator, R_estimator, v(k_start:k_step:k_end),...
         om(k_start:k_step:k_end), r(k_start:k_step:k_end,:),b(k_start:k_step:k_end,:), t(k_start:k_step:k_end),d,...
         fixed_landmarks, x_true(k_start:k_step:k_end), y_true(k_start:k_step:k_end), th_true(k_start:k_step:k_end),l);
    %toc


    % Landmark_estimates = X_estimate(end-1:end);
    % Landmark_estimates_Con = X_estimate_Con(end-1:end);
    % Landmark_estimates_Pen = X_estimate_Pen(end-1:end);
    %  X_estimate = X_estimate(1:end-2);
    % X_estimate_Con = X_estimate_Con(1:end-2);
    % X_estimate_Pen = X_estimate_Pen(1:end-2);
    %%
    error_calc;
    toc
    
    MC_errors(counter,:) = [rms_batch, rms_eq, rms_prob, rms_var_rem];
    disp('Done!')
end

% save('results_noiseless2.mat','MC_errors');
    %plotting;