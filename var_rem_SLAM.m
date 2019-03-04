function [X_op,X_L_op, rms_error_list, P_states, P_landmarks]  = var_rem_SLAM( X_0,P_0,X_L_0, Q, R, Y_v,Y_w, Y_r,Y_b, t, d,constr_L, x_true, y_true, th_true, l_true)
%GN_ESTIMATOR Summary of this function goes here
%% Some matrix dimensions
M = 17*2; %No of landmarks times dimension of landmark position in the state
dim = length(X_0); %No of trajectory state variables at timestep k
K = dim*(length(t)); %Dimension of the full trajectory
N = 2*nnz(Y_r); % Dimension of all available range and bearing measurements
T = t(2) - t(1);

%% For Plotting
rms_error_list = [];

%% Non-Linear Process Model
f_k = @(X, U) X + [T*cos(X(3)), 0;
                T*sin(X(3)), 0;
                0,        T]*U;
F_k = @(X_op, U) [1, 0, -T*U(1)*sin(X_op(3));
                  0, 1,  T*U(1)*cos(X_op(3));
                  0, 0, 1];
L_k = @(X_op) [T*cos(X_op(3)), 0;
               T*sin(X_op(3)), 0;
               0, T];
           



%% Non-linear Measurement Model

g_r_k = @(X_L, X) sqrt((X_L(2) - X(2)-d*sin(X(3))).^2 + (X_L(1) - X(1)-d*cos(X(3))).^2);
g_b_k = @(X_L, X) atan2(X_L(2) - X(2)-d*sin(X(3)), X_L(1) - X(1)-d*cos(X(3))) - X(3);



G1_r_mk = @(X_L, X) [X(1) + d*cos(X(3)) - X_L(1), X(2) + d*sin(X(3)) - X_L(2),...
          ((d*sin(X(3)))*(- X(1) - d*cos(X(3)) + X_L(1)) -(d*cos(X(3)))*(-X(2) - d*sin(X(3)) + X_L(2))) ]./g_r_k(X_L,X);

                                                                        
G1_b_mk = @(X_L, X) [X_L(2) - X(2) - d*sin(X(3)), -X_L(1) + X(1) + d*cos(X(3)),...
                      ((d*sin(X(3))*(-X(2) - d*sin(X(3)) + X_L(2))) + (d*cos(X(3))*(-X(1) - d*cos(X(3)) + X_L(1)))) - (g_r_k(X_L,X).^2)]./(g_r_k(X_L,X).^2);  

G1_mk = @(X_L,X) [G1_r_mk(X_L,X); G1_b_mk(X_L,X)];
    
G2_r_mk = @(X_L, X) [X_L(1) - X(1) - d*cos(X(3)), X_L(2) - X(2) - d*sin(X(3))]./g_r_k(X_L,X);
G2_b_mk = @(X_L, X) [-X_L(2) + X(2) + d*sin(X(3)), X_L(1) - X(1) - d*cos(X(3))]./g_r_k(X_L,X).^2;

G2_mk = @(X_L, X) [G2_r_mk(X_L,X); G2_b_mk(X_L,X)];


%% Initial States

X_op = randn(K,1);
X_op(1:dim) = X_0;
meas_no = 1;
for i=2*dim:dim:K
    X_op(i-(dim-1):i) = f_k(X_op(i-(2*dim - 1):i-dim), [Y_v(meas_no);Y_w(meas_no)]);
    X_op(i) = wrapToPi(X_op(i));
    meas_no = meas_no + 1;
end
% X_op(1:3:end) = x_true;
% X_op(2:3:end) = y_true;
% X_op(3:3:end) = th_true;
X_L_op = X_L_0;

%% Main iterative loop
 %%Allocating space for the various sparse matrices
    H_motion = spalloc(K,K, (K + dim*(K-dim)));
    H_obs = spalloc(N,K,dim*N);
    H_obs_2 = zeros(N,M);
    W_inv_motion = spalloc(K,K,dim*K);
    W_obs_inv = sparse(1:N,1:N,2*N);
    eps = 0.000001;
    
bar = waitbar(0,'Variable Removal SLAM iteration #1');

while(1)
 %% System Information Matrix
    
   
    %%Column by column create the H and W matrices corresponding to the process
    meas_no = 1;
    for k = 1:dim:K-dim
       %%Sparse matrix stuff

       %% Creating the H and W matrices
       temp = sparse([eye(dim);-F_k(X_op(k:k+dim-1),[Y_v(meas_no);Y_w(meas_no)])]);
       H_motion(k:k+2*dim-1,k:k+dim-1) = temp;
       temp = sparse(L_k(X_op(k:k+dim-1))*Q*L_k(X_op(k:k+dim-1)).') + eps*speye(dim);
       W_inv_motion(k:k+dim-1,k:k+dim-1) = inv(temp); 
       meas_no = meas_no + 1;
    end
    H_motion(K-(dim-1):K,K-(dim-1):K) = speye(dim);
    W_inv_motion(K-(dim-1):K,K-(dim-1):K) = sparse(inv(L_k(X_op(K-(dim-1):K))*Q*L_k(X_op(K-(dim-1):K)).' + eps*speye(dim)));
    W_inv_motion(1:dim,1:dim) = inv(P_0);
    
    
    %W_inv_motion = W_inv_motion + eps*speye(K);
    

    %%Column by Column creating the H and W matrices corresponding to
    %%measurements
    row = 1;
    col_s = 1;
    
    for j = 1:length(t) %Loop through every time step
         col_L = 1;
        for m = 1:M/2     %Loop through every measurement at a time step
            if Y_r(j,m) == 0 %checks if range measurement to landmark m is available
                col_L = col_L + 2;
                continue
            else
                H_obs(row:row+1,col_s:col_s + dim-1) = G1_mk(X_L_op(col_L:col_L+1), X_op(col_s:col_s+dim-1));
                H_obs_2(row:row+1,col_L:col_L + 1) = G2_mk(X_L_op(col_L:col_L+1), X_op(col_s:col_s+dim-1));
                W_obs_inv(row:row+1,row:row+1) = inv(R);
                row = row + 2;
                col_L = col_L + 2;
                
            end
        end
        col_s = col_s + dim;
        
    end
    
    
    H = [H_motion,zeros(K,M); H_obs, H_obs_2];
    W_inv = [W_inv_motion, sparse(K,N); sparse(N,K), W_obs_inv];
    
    System_matrix = (H.'*W_inv*H);
    
    %% System Information Vector

    e_op_X = zeros(K,1);
    e_op_X(1:dim) = X_0(1:dim) - X_op(1:dim);
%     e_op_X(3) = wrapToPi(e_op_X(3));
 

    meas_no = 1;
    for i=2*dim:dim:K
        predicted_X = f_k(X_op(i-(2*dim - 1):i-dim), [Y_v(meas_no);Y_w(meas_no)]);
        predicted_X(3) = wrapToPi(predicted_X(3));
%         e_op_X(i-(dim-1):i) = f_k(X_op(i-(2*dim - 1):i-dim), [Y_v(meas_no);Y_w(meas_no)]) - X_op(i-(dim-1):i);
        e_op_X(i-(dim-1):i) = predicted_X - X_op(i-(dim-1):i);
        if(abs(e_op_X(i))) > 5.0
            e_op_X(i) = predicted_X(3) + X_op(i);
        end
%         e_op_X(i) = wrapToPi(e_op_X(i));
        meas_no = meas_no + 1;
    end

    e_op_y = zeros(N,1);

    row = 1;
    for i = 1:length(t)
        for m = 1:M/2
            if Y_r(i,m) == 0
                continue
            else
                predicted_y_r = g_r_k(X_L_op(2*m -1:2*m),X_op(dim*i-2:dim*i));
                predicted_y_b =  wrapToPi(g_b_k(X_L_op(2*m -1:2*m),X_op(dim*i-2:dim*i)));
                e_op_y(row) = Y_r(i,m) - predicted_y_r;
                e_op_y(row+1) = Y_b(i,m) - predicted_y_b;
%                 e_op_y(row) = Y_r(i,m) - g_r_k(X_L_op(2*m -1:2*m),X_op(dim*i-2:dim*i));
%                 e_op_y(row+1) = wrapToPi(Y_b(i,m) - g_b_k(X_L_op(2*m -1:2*m),X_op(dim*i-2:dim*i)));
                row = row + 2;
            end
        end
    end
    
    e_op = [ e_op_X; e_op_y];
    System_vector = H.'*W_inv*e_op;

    %% Gauss-Newton Update
    
    System_matrix = System_matrix(4:end,4:end);
    System_vector = System_vector(4:end);
    
    %Removing the constrained landmarks
    ind_removed = [2*constr_L-1; 2*constr_L];
    ind_removed = ind_removed + K-dim;

    
    System_matrix = removerows(System_matrix,ind_removed);
    System_matrix = System_matrix.';
    System_matrix = removerows(System_matrix,ind_removed);
    System_matrix = System_matrix.';
    System_vector = removerows(System_vector,ind_removed);
    
    
    step = 1.0;

    del_X = System_matrix\System_vector;
    
    mag_del_X = norm(del_X);
%     del_X_list = [del_X_list; mag_del_X]; % For plotting
    
    X_op(4:end) = X_op(4:end) + step*del_X(1:end-(M-length(ind_removed)));
    X_op(3:3:end) = wrapToPi(X_op(3:3:end));
    
    del_X_L = zeros(M,1);
    counter_constr = 1;
    counter_meas = 1;
    for j=1:17
        if j == constr_L(counter_constr)
            del_X_L(2*j-1:2*j) = [0;0];
            counter_constr = counter_constr + 1;
        else            
            del_X_L(2*j-1:2*j) = del_X(K-dim + 2*counter_meas-1: K-dim + 2*counter_meas);
            counter_meas = counter_meas + 1;
        end
    end
    
    X_L_op(1:end) = X_L_op(1:end) + step*del_X_L;
    
    [x_var_rem, y_var_rem, l_var_rem] = align(X_op,X_L_op, x_true,y_true,l_true);
    
    rms_error = norm([x_true-x_var_rem; y_true - y_var_rem; l_true(:,1) - l_var_rem(:,1); l_true(:,2) - l_var_rem(:,2)]);
    rms_error_list = [rms_error_list; rms_error]; % For plotting
    
    if(mag_del_X<0.001 || length(rms_error_list)>100)
        waitbar(1,bar, 'Done!');
        pause(1);
        close(bar);
        %P_states = diag(inv(System_matrix(1:end-2,1:end-2)));
        %P_landmarks = diag(inv(System_matrix(end-1:end,end-1:end)));
        break
    end
    waitbar(length(rms_error_list)/200, bar, strcat('Variable Removal SLAM Iteration #',num2str(length(rms_error_list)))); 
end
    %P = speye(K,K)/System_matrix;
