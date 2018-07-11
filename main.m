%Batch LQG Estimation Exercise 2
%Author: Shatil Rahman

clear all
close all
clc
load dataset2.mat
%%

row = 1;
k_start = 1000;
k_end = 2000;

X_0 = [x_true(k_start); y_true(k_start); wrapToPi(th_true(k_start))]; 
P_0 = [0.01, 0.0, 0.0;
       0.0, 0.01, 0.0
       0.0, 0.0, 0.001];  

X_L_0 = 3*randn(2*length(l),1);
for i=1:length(l)
    X_L_0(row) = l(i,1) + 50;
    X_L_0(row+1) = l(i,2) + 50;
    row = row + 2;
end


Q_estimator = [v_var, 0;
               0, om_var];
R_estimator = [r_var, 0;
                0, b_var];
constrained_landmarks = [17,2;  2,7; 7,3; 17,7];
Y_const = zeros(length(constrained_landmarks),1);
for i = 1:length(constrained_landmarks)
   Y_const(i) = norm(l(constrained_landmarks(i,1),:) -  l(constrained_landmarks(i,2),:));
end

%%
tic
[q_batch,L_batch, del_Xs_batch] = GN_Estimator( X_0, P_0,X_L_0, Q_estimator, R_estimator, v(k_start:k_end),...
     om(k_start:k_end), r(k_start:k_end,:),b(k_start:k_end,:), t(k_start:k_end),d);
toc
tic
[q_constr,L_constr, del_Xs_constr] = ineq_SLAM2( X_0, P_0,X_L_0, Q_estimator, R_estimator, v(k_start:k_end),...
    om(k_start:k_end), r(k_start:k_end,:),b(k_start:k_end,:),Y_const, t(k_start:k_end),d,...
       constrained_landmarks);
toc
 
tic
[q_penalty,L_penalty, del_Xs_penalty] = constrained_SLAM( X_0, P_0,X_L_0, Q_estimator, R_estimator, v(k_start:k_end),...
    om(k_start:k_end), r(k_start:k_end,:),b(k_start:k_end,:),Y_const, t(k_start:k_end),d,...
       constrained_landmarks);
toc
% Landmark_estimates = X_estimate(end-1:end);
% Landmark_estimates_Con = X_estimate_Con(end-1:end);
% Landmark_estimates_Pen = X_estimate_Pen(end-1:end);
%  X_estimate = X_estimate(1:end-2);
% X_estimate_Con = X_estimate_Con(1:end-2);
% X_estimate_Pen = X_estimate_Pen(1:end-2);
%%
  
plotting;