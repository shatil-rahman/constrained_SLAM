%% Re-aligning the estimates for comparison

[x_batch, y_batch, l_batch] = align(q_batch,L_batch, x_true(k_start:k_step:k_end),y_true(k_start:k_step:k_end),l);
[x_prob, y_prob, l_prob] = align(q_prob,L_prob, x_true(k_start:k_step:k_end),y_true(k_start:k_step:k_end),l);
 [x_eq, y_eq, l_eq] = align(q_eq,L_eq, x_true(k_start:k_step:k_end),y_true(k_start:k_step:k_end),l);
[x_var_rem, y_var_rem, l_var_rem] = align(q_var_rem,L_var_rem, x_true(k_start:k_step:k_end),y_true(k_start:k_step:k_end),l);

%% RMS error calculations
rms_batch = norm([x_true(k_start:k_step:k_end)-x_batch; y_true(k_start:k_step:k_end) - y_batch; l(:,1) - l_batch(:,1); l(:,2) - l_batch(:,2)]);
rms_prob = norm([x_true(k_start:k_step:k_end)-x_prob; y_true(k_start:k_step:k_end) - y_prob; l(:,1) - l_prob(:,1); l(:,2) - l_prob(:,2)]);
rms_eq = norm([x_true(k_start:k_step:k_end)-x_eq; y_true(k_start:k_step:k_end) - y_eq; l(:,1) - l_eq(:,1); l(:,2) - l_eq(:,2)]);
rms_var_rem = norm([x_true(k_start:k_step:k_end)-x_var_rem; y_true(k_start:k_step:k_end) - y_var_rem; l(:,1) - l_var_rem(:,1); l(:,2) - l_var_rem(:,2)]);

