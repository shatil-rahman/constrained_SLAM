


%% Re-aligning the estimates for comparison

[x_batch, y_batch, l_batch] = align(q_batch,L_batch, x_true(k_start:k_step:k_end),y_true(k_start:k_step:k_end),l);
[x_constr, y_constr, l_constr] = align(q_constr,L_constr, x_true(k_start:k_step:k_end),y_true(k_start:k_step:k_end),l);
[x_penalty, y_penalty, l_penalty] = align(q_penalty,L_penalty, x_true(k_start:k_step:k_end),y_true(k_start:k_step:k_end),l);

%% RMS error calculations
rms_batch = norm([x_true(k_start:k_step:k_end)-x_batch; y_true(k_start:k_step:k_end) - y_batch; l(:,1) - l_batch(:,1); l(:,2) - l_batch(:,2)]);
rms_constr = norm([x_true(k_start:k_step:k_end)-x_constr; y_true(k_start:k_step:k_end) - y_constr; l(:,1) - l_constr(:,1); l(:,2) - l_constr(:,2)]);
rms_penalty = norm([x_true(k_start:k_step:k_end)-x_penalty; y_true(k_start:k_step:k_end) - y_penalty; l(:,1) - l_penalty(:,1); l(:,2) - l_penalty(:,2)]);
% Plotting

figure
plot(x_true(k_start:k_step:k_end),y_true(k_start:k_step:k_end),'-')
title('SLAM','Interpreter','laTex')
xlabel('x position (m)','Interpreter','laTex')
ylabel('y position, (m)','Interpreter','laTex')
grid on
hold on
plot(x_batch,y_batch,'g-')
plot(x_constr,y_constr,'m-')
% plot(x_penalty,y_penalty,'r-')
plot(l(:,1),l(:,2),'bo')
plot(l_batch(:,1),l_batch(:,2),'go')
plot(l_constr(:,1),l_constr(:,2),'mo')
% plot(l_penalty(:,1),l_penalty(:,2),'ro')

% plot(t,X_estimate_Con(1:2:end),'-')
% plot(t,X_estimate_Pen(1:2:end),'-')
legend({'Groundtruth','Unconstrained Estimate','Stochastically Constrained Estimate'},'Interpreter','laTex')%, 'True Landmark Position', 'Batch Landmark Position', 'Penalty Landmark Position')
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
% plot(del_Xs_penalty)
% plot(del_Xs_Pen)
title('RMS Error vs Number of Iterations','Interpreter','laTex')
xlabel('No of iterations','Interpreter','laTex')
ylabel('RMS Error','Interpreter','laTex')
legend({'Unconstrained Estimate','Stochastically Constrained Estimate'},'Interpreter','laTex')

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