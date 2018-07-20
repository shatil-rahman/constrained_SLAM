function [Cqp, std_dev] = DCM_estimation( d_p, d_q,odom,odom_var,T,w_p, w_q,pairs_p, pairs_q )
%DCM_ESTIMATION Finds Cqp such that d_q = Cqp*dp

    %% If one of the poses has only seen 1 landmark, use odometry only
    
    if isempty(d_p) || isempty(d_q)
       Cqp = [cos(T*odom), sin(T*odom), 0;
              -sin(T*odom), cos(T*odom),0;
              0, 0, 1];
       std_dev = odom_var;
       return
    end


    %% Picking out the common relative feature measurements
    [l_pq, ip, iq] = intersect(pairs_p, pairs_q, 'rows');
    eps = 0.000001;
    
    d_p_com = zeros(3*length(ip),1);
    d_q_com = zeros(3*length(ip),1);
    d_p_com(1:3:end) = d_p(ip,1);
    d_p_com(2:3:end) = d_p(ip,2);
    d_q_com(1:3:end) = d_q(iq,1);
    d_q_com(2:3:end) = d_q(iq,2);
    
    
    
     w_bar = zeros(length(ip)*3,length(ip)*3);
     
     for i = 1:length(ip)
         w_p_i = [w_p(2*i-1:2*i,2*i-1:2*i), zeros(2,1); zeros(1,2),eps];
         w_q_i = [w_q(2*i-1:2*i,2*i-1:2*i), zeros(2,1); zeros(1,2),eps];
         w_bar(3*i-2:3*i,3*i-2:3*i) = w_q_i + w_p_i - ((w_p_i - w_q_i).'/(w_q_i + w_p_i))*((w_p_i - w_q_i));
     end
%     w_p_com = zeros(length(ip)*3,length(ip)*3);
%     w_p_com(3.*(1:length(ip))-2,3.*(1:length(ip))-2) = w_p(2.*ip-1,2.*ip-1);
%     w_p_com(3.*(1:length(ip))-1,3.*(1:length(ip))-1) = w_p(2.*ip,2.*ip);
%     
%     w_q_com = zeros(length(ip)*3,length(ip)*3);
%     w_q_com(3.*(1:length(iq))-2,3.*(1:length(iq))-2) = w_q(2.*iq-1,2.*iq-1);
%     w_q_com(3.*(1:length(iq))-1,3.*(1:length(iq))-1) = w_q(2.*iq,2.*iq);
%     
%     w_q_com = w_q_com + eps*eye(length(w_q_com)); %SPD hack
%     w_p_com = w_p_com + eps*eye(length(w_p_com)); %SPD hack
    
    
    %% 
%     
%     [V_p, D_p ]= eig(w_p_com);
%     [V_q, D_q ]= eig(w_q_com);
%     
%     d_p_com = V_p.'*d_p_com;
%     d_q_com = V_q.'*d_q_com;
%     
%     w_p_com = D_p;
%     w_q_com = D_q;
    
    B = d_q_com - d_p_com;
    
    A = zeros(length(d_p_com),3);
    for i=1:length(ip)
        temp = d_q_com(3*i-2:3*i) + d_p_com(3*i-2:3*i);
        A(3*i-2:3*i,:) = [0, -temp(3), temp(2);
                          temp(3), 0, -temp(1);
                          -temp(2), temp(1), 0];
                      
    end
    
%     W = w_q_com + w_p_com - ((w_p_com - w_q_com).'/(w_q_com + w_p_com))*((w_p_com - w_q_com));
%     W = W + eps*eye(length(W));
%     W = eye(length(W));
    
    w_odom = [eps, 0, 0;
                0, eps, 0;
                0, 0, odom_var];
    A_odom = eye(3);
    B_odom = [0;0;tan(T*odom/2)];
    
    A = [A_odom; A];
    B = [B_odom; B];
    w_bar = [w_odom, zeros(3,length(w_bar));
             zeros(length(w_bar),3), w_bar;];

    A_k = ((A.'/w_bar)*A);
    B_k = (A.'/w_bar)*B;

    A_k = A_k + eye(length(A_k)).*eps;
    
%     A_k = [A_k; A_odom];
%     B_k = [B_k; B_odom];
    p = A_k\B_k;
    
    temp = [0, -p(3), p(2);
                p(3), 0, -p(1);
                -p(2), p(1), 0];
    Cqp = (eye(3) + temp)\((eye(3)- temp));
    covariance = inv(A_k);
    sigma_3 = covariance(3,3);
    std_dev = ((2/(1+p(3)^2))^2)*sigma_3;

end

