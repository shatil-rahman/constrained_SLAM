function [ X_est_tr, Y_est_tr, L_est_tr ] = align(X_estimate, L_estimate, x_true, y_true, L_true)
%ALIGN Uses Horn's method to align the SLAM estimate with true global frame
%   Detailed explanation goes here
k = length(x_true);

A = [X_estimate(1:3:end).', L_estimate(1:2:end).';X_estimate(2:3:end).',L_estimate(2:2:end).'];
B = [x_true.',L_true(:,1).'; y_true.',L_true(:,2).'];

[reg_params, transformed_estimates] = absor(A,B);
X_est_tr = transformed_estimates(1,1:k);
Y_est_tr = transformed_estimates(2,1:k);
L_est_tr = [transformed_estimates(1,k+1:end).',transformed_estimates(2,k+1:end).'];

end

