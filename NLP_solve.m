%% Script to call the constrained optimizer and related functions

global H;
global W_inv;
load dataset2.mat;


%%Initialize states and landmarks
T = t(2) - t(1);
X_0 = [x_true(1); y_true(1); wrapToPi(th_true(1))];

f_k = @(X, U) X + [T*cos(X(3)), 0;
                T*sin(X(3)), 0;
                0,        T]*U; %process model, integrating to initialize states
dim = 3;
k_end = 2000;

K = dim*length(t(1:k_end));            
X_op = ones(K,1)*10;
X_op(1:dim) = X_0;
meas_no = 1;
for i=2*dim:dim:K
    X_op(i-(dim-1):i) = f_k(X_op(i-(2*dim - 1):i-dim), [v(meas_no);om(meas_no)]);
    X_op(i) = wrapToPi(X_op(i));
    meas_no = meas_no + 1;
end            


