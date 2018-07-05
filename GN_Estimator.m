function [X_op,X_L_op, del_X_list, P_states, P_landmarks]  = GN_Estimator( X_0,P_0,X_L_0, Q, R, Y_v,Y_w, Y_r,Y_b, t, d)
%GN_ESTIMATOR Summary of this function goes here
%% Some matrix dimensions
M = 17*2; %No of landmarks times dimension of landmark position in the state
dim = length(X_0); %No of trajectory state variables at timestep k
K = dim*(length(t)); %Dimension of the full trajectory
N = 2*nnz(Y_r); % Dimension of all available range and bearing measurements
T = t(2) - t(1);

%% For Plotting
del_X_list = [];

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
g_b_k = @(X_L, X) wrapToPi(atan2(X_L(2) - X(2)-d*sin(X(3)), X_L(1) - X(1)-d*cos(X(3))) - X(3));



G1_r_mk = @(X_L, X) [X(1) + d*cos(X(3)) - X_L(1), X(2) + d*sin(X(3)) - X_L(2),...
          ((d*sin(X(3)))*(- X(1) - d*cos(X(3)) + X_L(1)) -(d*cos(X(3)))*(-X(2) - d*sin(X(3)) + X_L(2))) ]./g_r_k(X_L,X);

                                                                        
G1_b_mk = @(X_L, X) [X_L(2) - X(2) - d*sin(X(3)), -X_L(1) + X(1) + d*cos(X(3)),...
                      ((d*sin(X(3))*(-X(2) - d*sin(X(3)) + X_L(2))) + (d*cos(X(3))*(-X(1) - d*cos(X(3)) + X_L(1)))) - (g_r_k(X_L,X).^2)]./(g_r_k(X_L,X).^2);  

G1_mk = @(X_L,X) [G1_r_mk(X_L,X); G1_b_mk(X_L,X)];
    
G2_r_mk = @(X_L, X) [X_L(1) - X(1) - d*cos(X(3)), X_L(2) - X(2) - d*sin(X(3))]./g_r_k(X_L,X);
G2_b_mk = @(X_L, X) [-X_L(2) + X(2) + d*sin(X(3)), X_L(1) - X(1) - d*cos(X(3))]./g_r_k(X_L,X).^2;

G2_mk = @(X_L, X) [G2_r_mk(X_L,X); G2_b_mk(X_L,X)];


%% Initial States

X_op = ones(K,1)*10;
X_op(1:dim) = X_0;
meas_no = 1;
for i=2*dim:dim:K
    X_op(i-(dim-1):i) = f_k(X_op(i-(2*dim - 1):i-dim), [Y_v(meas_no);Y_w(meas_no)]);
    X_op(i) = wrapToPi(X_op(i));
    meas_no = meas_no + 1;
end

X_L_op = X_L_0;


%% Main iterative loop
 %%Allocating space for the various sparse matrices
    H_motion = spalloc(K,K, (K + dim*(K-dim)));
    H_obs = spalloc(N,K,dim*N);
    H_obs_2 = zeros(N,M);
    W_inv_motion = spalloc(K,K,dim*K);
    W_obs_inv = sparse(1:N,1:N,2*N);
    eps = 0.000001;
    
bar = waitbar(0,'Batch SLAM iteration #1');

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
    e_op_X(3) = wrapToPi(e_op_X(3));
 

    meas_no = 1;
    for i=2*dim:dim:K
%         predicted_X = f_k(X_op(i-(2*dim - 1):i-dim), [Y_v(meas_no);Y_w(meas_no)]);
%         predicted_X(3) = wrapToPi(predicted_X(3));
        e_op_X(i-(dim-1):i) = f_k(X_op(i-(2*dim - 1):i-dim), [Y_v(meas_no);Y_w(meas_no)]) - X_op(i-(dim-1):i);
        e_op_X(i) = wrapToPi(e_op_X(i));
        meas_no = meas_no + 1;
    end

    e_op_y = zeros(N,1);

    row = 1;
    for i = 1:length(t)
        for m = 1:M/2
            if Y_r(i,m) == 0
                continue
            else
                e_op_y(row) = Y_r(i,m) - g_r_k(X_L_op(2*m -1:2*m),X_op(dim*i-2:dim*i));
                e_op_y(row+1) = wrapToPi(Y_b(i,m) - g_b_k(X_L_op(2*m -1:2*m),X_op(dim*i-2:dim*i)));
                row = row + 2;
            end
        end
    end
    
    e_op = [ e_op_X; e_op_y];
    System_vector = H.'*W_inv*e_op;

    %% Gauss-Newton Update
    
    System_matrix = System_matrix(4:end,4:end);
    System_vector = System_vector(4:end);

    del_X = System_matrix\System_vector;
    
    mag_del_X = norm(del_X);
    del_X_list = [del_X_list; mag_del_X]; % For plotting
    
    X_op(4:end) = X_op(4:end) + del_X(1:end-M);
    X_op(3:3:end) = wrapToPi(X_op(3:3:end));
    X_L_op(1:end) = X_L_op(1:end) + del_X(end-M+1:end);
    
    if(mag_del_X<0.001 || length(del_X_list)>20)
        waitbar(1,bar, 'Done!');
        pause(1);
        close(bar);
        %P_states = diag(inv(System_matrix(1:end-2,1:end-2)));
        %P_landmarks = diag(inv(System_matrix(end-1:end,end-1:end)));
        break
    end
    waitbar(length(del_X_list)/20, bar, strcat('Batch SLAM Iteration #',num2str(length(del_X_list)))); 
end
    %P = speye(K,K)/System_matrix;
