function [Cqp] = DCM_estimation( d_p, d_q, w_p, w_q, pairs_p, pairs_q )
%DCM_ESTIMATION Finds Cqp such that d_q = Cqp*dp

    [l_pq, ip, iq] = intersect(pairs_p, pairs_q, 'rows');
    eps = 0.000001;
    
    d_p_com = zeros(3*length(l_pq),1);
    d_q_com = zeros(3*length(l_pq),1);
    d_p_com(1:3:end) = d_p(ip,1);
    d_p_com(2:3:end) = d_p(ip,2);
    d_q_com(1:3:end) = d_q(iq,1);
    d_q_com(2:3:end) = d_q(iq,2);
    
    w_p_com = zeros(length(l_pq)*3,length(l_pq)*3);
    w_p_com(3.*(1:length(ip))-2,3.*(1:length(ip))-2) = w_p(2.*ip-1,2.*ip-1);
    w_p_com(3.*(1:length(ip))-1,3.*(1:length(ip))-1) = w_p(2.*ip,2.*ip);
    
    w_q_com = zeros(length(l_pq)*3,length(l_pq)*3);
    w_q_com(3.*(1:length(iq))-2,3.*(1:length(iq))-2) = w_q(2.*iq-1,2.*iq-1);
    w_q_com(3.*(1:length(iq))-1,3.*(1:length(iq))-1) = w_q(2.*iq,2.*iq);
    
    w_q_com = w_q_com + eps*eye(length(w_q_com));
    w_p_com = w_p_com + eps*eye(length(w_p_com));
    
    B = d_q_com - d_p_com;
    
    A = zeros(length(d_p_com),3);
    for i=1:length(l_pq)
        temp = d_q_com(3*i-2:3*i) + d_p_com(3*i-2:3*i);
        A(3*i-2:3*i,:) = [0, -temp(3), temp(2);
                          temp(3), 0, -temp(1);
                          -temp(2), temp(1), 0];
                      
    end
    
    W = w_q_com + w_p_com - ((w_p_com - w_q_com).'/(w_q_com + w_p_com))*((w_p_com - w_q_com));
    
    A_k = (A.'/W*A);
    B_k = A.'/W*B;
    
    p = A_k\B_k;

end

