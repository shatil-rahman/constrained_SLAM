function [X_op,X_L_op, del_X_list, P_states, P_landmarks]  = ineq_SLAM( X_0,P_0,X_L_0, Q, R, Y_v,Y_w, Y_r,Y_b,Y_const, t, d,constr_L)
%GN_ESTIMATOR Summary of this function goes here
%% Some matrix dimensions
M = 17*2; %No of landmarks times dimension of landmark position in the state
dim = length(X_0); %No of trajectory state variables at timestep k
K = dim*(length(t)); %Dimension of the full trajectory
N = 2*nnz(Y_r); % Dimension of all available range and bearing measurements
T = t(2) - t(1);
C = length(Y_const);

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
%g_b_k = @(X_L, X) atan2(X_L(2) - X(2)-d*sin(X(3)), X_L(1) - X(1)-d*cos(X(3))) - X(3);



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
meas_n = 1;
for i=2*dim:dim:K
    X_op(i-(dim-1):i) = f_k(X_op(i-(2*dim - 1):i-dim), [Y_v(meas_n);Y_w(meas_n)]);
    X_op(i) = wrapToPi(X_op(i));
    meas_n = meas_n + 1;
end

X_L_op = X_L_0;


 %%Allocating space for the various sparse matrices
    H_motion = spalloc(K,K, (K + dim*(K-dim)));
    H_obs = spalloc(N,K,dim*N);
    H_obs_2 = zeros(N,M);
    W_inv_motion = spalloc(K,K,dim*K);
    W_obs_inv = sparse(1:N,1:N,2*N);
    W_inv = [];
    eps = 0.000001;
    eps_bound = 1.0;
    X_last = [];
    gradJ = [];
    HessJ = [];
    

x0 = [X_op; X_L_op];
options = optimoptions('fmincon','Algorithm','interior-point','SpecifyObjectiveGradient',true,...
                        'SpecifyConstraintGradient',true,'HessianFcn',@hessianfcn,'PlotFcn','optimplotfval',...
                        'Display','iter-detailed');
                    
[X, fval, exitflag, output] = fmincon(@objectiveFunc, x0,[],[],[],[],[],[],@constraints,options);
                    
    %% Function definitions
    
    function [System_vector, System_matrix] = updateGradHess(X_op, X_L_op)
     
        X_op(3:3:end) = wrapToPi(X_op(3:3:end));
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
        System_matrix = System_matrix(4:end,4:end);

        %% System Information Vector

        e_op_X = zeros(K,1);
        e_op_X(1:dim) = X_0(1:dim) - X_op(1:dim);
        e_op_X(3) = wrapToPi(e_op_X(3));


        meas_no = 1;
        for j=2*dim:dim:K

            e_op_X(j-(dim-1):j) = f_k(X_op(j-(2*dim - 1):j-dim), [Y_v(meas_no);Y_w(meas_no)]) - X_op(j-(dim-1):j);
            e_op_X(j) = wrapToPi(e_op_X(j));
            meas_no = meas_no + 1;
        end

        e_op_y = zeros(N,1);

        row = 1;
        for j = 1:length(t)
            for m = 1:M/2
                if Y_r(j,m) == 0
                    continue
                else
                    e_op_y(row) = Y_r(j,m) - g_r_k(X_L_op(2*m -1:2*m),X_op(dim*j-2:dim*j));
                    e_op_y(row+1) = wrapToPi(Y_b(j,m) - g_b_k(X_L_op(2*m -1:2*m),X_op(dim*j-2:dim*j)));
                    row = row + 2;
                end
            end
        end

        e_op = [ e_op_X; e_op_y];
        System_vector = -H.'*W_inv*e_op;
        System_vector = System_vector(4:end);
    end

    function [J, grad] = objectiveFunc(X)
     
        X_state = X(1:end-M);
        X_state(3:3:end) = wrapToPi(X_state(3:3:end));
        X_L = X(end-M+1:end);
        e_op_X = zeros(K,1);
        e_op_X(1:dim) = X_0(1:dim) - X_state(1:dim);
        e_op_X(3) = wrapToPi(e_op_X(3));
        
        if ~isequal(X,X_last)
            [gradJ, HessJ] = updateGradHess(X_state, X_L);
            X_last = X;
        end
        grad = gradJ;

        meas_no = 1;
        for j=2*dim:dim:K

            e_op_X(j-(dim-1):j) = f_k(X_state(j-(2*dim - 1):j-dim), [Y_v(meas_no);Y_w(meas_no)]) - X_state(j-(dim-1):j);
            e_op_X(j) = wrapToPi(e_op_X(j));
            meas_no = meas_no + 1;
        end

        e_op_y = zeros(N,1);

        row = 1;
        for j = 1:length(t)
            for m = 1:M/2
                if Y_r(j,m) == 0
                    continue
                else
                    e_op_y(row) = Y_r(j,m) - g_r_k(X_L(2*m -1:2*m),X_state(dim*j-2:dim*j));
                    e_op_y(row+1) = wrapToPi(Y_b(j,m)) - wrapToPi(g_b_k(X_L(2*m -1:2*m),X_state(dim*j-2:dim*j)));
                    %e_op_y(row+1) = Y_b(j,m) - g_b_k(X_L(2*m -1:2*m),X_state(dim*j-2:dim*j));
                    row = row + 2;
                end
            end
        end

        e_op = [ e_op_X; e_op_y];
        e_op = e_op(4:end);
        
        J = 0.5* e_op.'*W_inv(4:end,4:end)*e_op;
        
        
    end

    function [c, ceq, gc, geq] = constraints(X)
        dist = @(L1,L2) sqrt((L2(1)-L1(1)).^2 + (L2(2)-L1(2)).^2);
        
        X_state = X(1:end-M);
        X_L = X(end-M+1:end);
        
        c_no = 1;
        c = zeros(C*2,1);
        J_c = spalloc(C*2,length(X),4*C);
        for j = 1:C
            m_1= constr_L(j,1);
            m_2= constr_L(j,2);
            L1 = X_L(2*m_1-1:2*m_1);
            L2 = X_L(2*m_2-1:2*m_2);
            
            c(c_no) = dist(L1, L2) - Y_const(j) - eps_bound;
            c(c_no + 1) = -dist(L1, L2) + Y_const(j) - eps_bound;
            
            J_c(c_no,K + 2*m_1-1:K+2*m_1) = [L1(1) - L2(1), L1(2) - L2(2)]./dist(L1,L2);
            J_c(c_no,K + 2*m_2-1:K+2*m_2) = [L2(1) - L1(1), L2(2) - L1(2)]./dist(L1,L2);
            J_c(c_no+1,K + 2*m_1-1:K+2*m_1) = -J_c(c_no,K + 2*m_1-1:K+2*m_1);
            J_c(c_no+1,K + 2*m_2-1:K+2*m_2) = -J_c(c_no,K + 2*m_2-1:K+2*m_2); 
            
            gc = J_c.';
            
            c_no = c_no + 2;
        end
        
        ceq = [];
        geq = [];
    end



    function hessian = hessianfcn(X,lambda)
        dist = @(L1,L2) sqrt((L2(1)-L1(1)).^2 + (L2(2)-L1(2)).^2);
        
        
        X_state = X(1:end-M);
        X_state(3:3:end) = wrapToPi(X_state(3:3:end));
        X_L = X(end-M+1:end);
        
        if ~isequal(X,X_last)
            [gradJ, HessJ] = updateGradHess(X_state, X_L);
            X_last = X;
        end
        hessian = HessJ;
        
        
        c_no = 1;

        for j=1:C
            m_1= constr_L(j,1);
            m_2= constr_L(j,2);
            L1 = X_L(2*m_1-1:2*m_1);
            L2 = X_L(2*m_2-1:2*m_2);
            
            dx1_dx1 = 1/dist(L1,L2) - (L1(1)-L2(1)).^2/(dist(L1,L2)^.3);
            dy1_dy1 = 1/dist(L1,L2) - (L1(2)-L2(2)).^2/(dist(L1,L2)^.3);
            dx1_dy1 = -((L1(1)-L2(1))*(L1(2)-L2(2)))/(dist(L1,L2).^3);
             
            
            %dx1 row
            hessian(K+2*m_1-1,K+2*m_1-1) = hessian(K+2*m_1-1,K+2*m_1-1) + lambda.ineqnonlin(c_no).*dx1_dx1 + lambda.ineqnonlin(c_no+1).*(-dx1_dx1);
            hessian(K+2*m_1-1,K+2*m_1)   = hessian(K+2*m_1-1,K+2*m_1) + lambda.ineqnonlin(c_no).*dx1_dy1+ lambda.ineqnonlin(c_no+1).*(-dx1_dy1);
            hessian(K+2*m_1-1,K+2*m_2-1) = hessian(K+2*m_1-1,K+2*m_2-1) + lambda.ineqnonlin(c_no).*(-dx1_dx1) + lambda.ineqnonlin(c_no+1).*dx1_dx1;
            hessian(K+2*m_1-1,K+2*m_2)   = hessian(K+2*m_1-1,K+2*m_2) + lambda.ineqnonlin(c_no).*(-dx1_dy1) + lambda.ineqnonlin(c_no+1).*dx1_dy1;
            
            %dy1 row
            hessian(K+2*m_1,K+2*m_1-1) = hessian(K+2*m_1,K+2*m_1-1) + lambda.ineqnonlin(c_no).*dx1_dy1 + lambda.ineqnonlin(c_no+1).*(-dx1_dy1);
            hessian(K+2*m_1,K+2*m_1)   = hessian(K+2*m_1,K+2*m_1) + lambda.ineqnonlin(c_no).*dy1_dy1+ lambda.ineqnonlin(c_no+1).*(-dy1_dy1);
            hessian(K+2*m_1,K+2*m_2-1) = hessian(K+2*m_1,K+2*m_2-1) + lambda.ineqnonlin(c_no).*(-dx1_dy1) + lambda.ineqnonlin(c_no+1).*dx1_dy1;
            hessian(K+2*m_1,K+2*m_2)   = hessian(K+2*m_1,K+2*m_2) + lambda.ineqnonlin(c_no).*(-dy1_dy1) + lambda.ineqnonlin(c_no+1).*dy1_dy1;
            
            %dx2 row
            hessian(K+2*m_2-1,K+2*m_1-1) = hessian(K+2*m_2-1,K+2*m_1-1) + -1*(lambda.ineqnonlin(c_no).*dx1_dx1 + lambda.ineqnonlin(c_no+1).*(-dx1_dx1));
            hessian(K+2*m_2-1,K+2*m_1) = hessian(K+2*m_2-1,K+2*m_1) + -1*(lambda.ineqnonlin(c_no).*dx1_dy1+ lambda.ineqnonlin(c_no+1).*(-dx1_dy1));
            hessian(K+2*m_2-1,K+2*m_2-1) = hessian(K+2*m_2-1,K+2*m_2-1) + -1*(lambda.ineqnonlin(c_no).*(-dx1_dx1) + lambda.ineqnonlin(c_no+1).*dx1_dx1);
            hessian(K+2*m_2-1,K+2*m_2) = hessian(K+2*m_2-1,K+2*m_2) + -1*(lambda.ineqnonlin(c_no).*(-dx1_dy1) + lambda.ineqnonlin(c_no+1).*dx1_dy1);
            
            %dy2 row
            hessian(K+2*m_2,K+2*m_1-1) = hessian(K+2*m_1,K+2*m_1-1) + -1*(lambda.ineqnonlin(c_no).*dx1_dy1 + lambda.ineqnonlin(c_no+1).*(-dx1_dy1));
            hessian(K+2*m_2,K+2*m_1)   = hessian(K+2*m_1,K+2*m_1) + -1*(lambda.ineqnonlin(c_no).*dy1_dy1+ lambda.ineqnonlin(c_no+1).*(-dy1_dy1));
            hessian(K+2*m_2,K+2*m_2-1) = hessian(K+2*m_1,K+2*m_2-1) + -1*(lambda.ineqnonlin(c_no).*(-dx1_dy1) + lambda.ineqnonlin(c_no+1).*dx1_dy1);
            hessian(K+2*m_2,K+2*m_2)   = hessian(K+2*m_1,K+2*m_2) + -1*(lambda.ineqnonlin(c_no).*(-dy1_dy1) + lambda.ineqnonlin(c_no+1).*dy1_dy1);
            
            c_no = c_no + 2;
        end
        
    end
end
    %P = speye(K,K)/System_matrix;

