%% Plotting
% pos_error = X_estimate(1:2:end) - X(1,:).';
% vel_error = X_estimate(2:2:end) - X(2,:).';
% pos_error_con = X_estimate_Con(1:2:end) - X(1,:).';
% vel_error_con = X_estimate_Con(2:2:end) - X(2,:).';
% pos_error_pen = X_estimate_Pen(1:2:end) - X(1,:).';
% vel_error_pen = X_estimate_Pen(2:2:end) - X(2,:).';
% three_sig_bounds = 3.*(sqrt(P_x));
% three_sig_bounds_con = 3.*(sqrt(P_x_Con));
% three_sig_bounds_pen = 3.*(sqrt(P_x_Pen));

%% Re-aligning the estimates for comparison

[x_batch, y_batch, l_batch] = align(q_batch,L_batch, x_true(k_start:k_end),y_true(k_start:k_end),l);
[x_constr, y_constr, l_constr] = align(q_constr,L_constr, x_true(k_start:k_end),y_true(k_start:k_end),l);
[x_penalty, y_penalty, l_penalty] = align(q_penalty,L_penalty, x_true(k_start:k_end),y_true(k_start:k_end),l);


%% Plotting

figure
plot(x_true(k_start:k_end),y_true(k_start:k_end),'-')
title('SLAM')
xlabel('x position (m)')
ylabel('y position, (m)')
grid on
hold on
plot(x_batch,y_batch,'g-')
plot(x_constr,y_constr,'m-')
plot(x_penalty,y_penalty,'r-')
plot(l(:,1),l(:,2),'bo')
plot(l_batch(:,1),l_batch(:,2),'go')
plot(l_constr(:,1),l_constr(:,2),'mo')
plot(l_penalty(:,1),l_penalty(:,2),'ro')

% plot(t,X_estimate_Con(1:2:end),'-')
% plot(t,X_estimate_Pen(1:2:end),'-')
legend('True Position','Batch Estimate','Constrained Estimate','Penalty Estimate')%, 'True Landmark Position', 'Batch Landmark Position', 'Penalty Landmark Position')
% subplot(2,1,2)
% plot(t,X(2,:),'-')
% title('Velocity')
% xlabel('Time (s)')
% ylabel('Velocity, (m/s)')
% grid on
% hold on
% plot(t,X_estimate(2:2:end),'-')
% plot(t,X_estimate_Con(2:2:end),'-')
% plot(t,X_estimate_Pen(2:2:end),'-')
% legend('True Velocity', 'Batch Estimate','Constrained Batch Estimate','Penalty Batch Estimate')

%% del_X vs time

figure
plot(del_Xs_batch)
hold on
plot(del_Xs_constr)
plot(del_Xs_penalty)
% plot(del_Xs_Pen)
title('\deltaX vs Number of Iterations')
xlabel('No of iterations')
ylabel('del_X')
legend('Regular Batch','Constrained Batch','Penalty Batch')

% %% Error plots
% 
% figure
% subplot(3,1,1)
% plot(t,pos_error_con,'-')
% hold on
% plot(t,three_sig_bounds_con(1:2:end),'g-')
% hold on
% plot(t,-1.*three_sig_bounds_con(1:2:end),'g-')
% grid on
% title(' Position Error')
% legend('Position Error', '3 sigma bound')
% 
% subplot(3,1,2)
% plot(t,pos_error_pen,'-')
% hold on
% plot(t,three_sig_bounds_pen(2:2:end),'g-')
% hold on
% plot(t,-1.*three_sig_bounds_pen(2:2:end),'g-')
% grid on
% title('Position Error - Penalty SLAM')
% legend('Position Error', '3 sigma bound')
% 
% subplot(3,1,3)
% plot(t,pos_error,'-')
% hold on
% plot(t,three_sig_bounds(2:2:end),'g-')
% hold on
% plot(t,-1.*three_sig_bounds(2:2:end),'g-')
% grid on
% title('Position Error - Batch SLAM')
% legend('Position Error', '3 sigma bound')