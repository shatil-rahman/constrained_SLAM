



%% Plotting

figure
plot(x_true(k_start:k_step:k_end),y_true(k_start:k_step:k_end),'--')
title('SLAM','Interpreter','laTex')
xlabel('x position (m)','Interpreter','laTex')
ylabel('y position, (m)','Interpreter','laTex')
grid on
hold on
plot(x_batch,y_batch,'g-')
plot(x_prob,y_prob,'r-')
plot(x_eq,y_eq,'m-')
plot(x_var_rem,y_var_rem,'k-')
plot(l(:,1),l(:,2),'bo')
plot(l_batch(:,1),l_batch(:,2),'go')
plot(l_prob(:,1),l_prob(:,2),'ro')
plot(l_eq(:,1),l_eq(:,2),'mo')
plot(l_var_rem(:,1),l_var_rem(:,2),'ko')

% plot(t,X_estimate_Con(1:2:end),'-')
% plot(t,X_estimate_Pen(1:2:end),'-')
legend({'Groundtruth','No constraints','Relative Probabilistic Constraints',...
    'Relative Equality Constraints','Variable Removal Method'},'Interpreter','laTex')%, 'True Landmark Position', 'Batch Landmark Position', 'Penalty Landmark Position')
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
plot(del_Xs_batch,'g-')
hold on
plot(del_Xs_prob,'r-')
plot(del_Xs_eq,'m-')
plot(del_Xs_var_rem,'k-')

title('RMS Error vs Number of Iterations','Interpreter','laTex')
xlabel('No of iterations','Interpreter','laTex')
ylabel('RMS Error','Interpreter','laTex')
legend({'No constraints','Relative Probabilistic Constraints',...
    'Relative Equality Constraints','Variable Removal Method'},'Interpreter','laTex')

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