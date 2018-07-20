function [ rel_disp, land_pairs, covariance ] = inv_meas(range, bearing, R_range,R_bearing)
%UNTITLED Computes the relative displacements between landmarks at a given
%for measurements at a given timestep
    if nnz(range)<2
        nnz(range)
        rel_disp = [];
        land_pairs = [];
        covariance = [];
        return
    end
        
    n_rel = nchoosek(nnz(range),2);
    land_pairs = zeros(n_rel,2);
    rel_disp = zeros(n_rel,2);
    Jac = zeros(n_rel*2, 17*2);
    meas_no = 1;
    eps = 0.00000001;

    for i = 1:length(range) - 1
        if range(i) == 0
            continue        
        else
            for j = i+1:length(range)
               if range(j) == 0
                   continue

               else
                   land_pairs(meas_no,:) =  [i,j];
                   d_x = range(j)*cos(bearing(j)) - range(i)*cos(bearing(i));
                   d_y = range(j)*sin(bearing(j)) - range(i)*sin(bearing(i));
                   rel_disp(meas_no,:) = [d_x, d_y];
                   
                   Jac(2*meas_no -1:2*meas_no,2*j-1:2*j) = [cos(bearing(j)), -range(j)*sin(bearing(j));...
                                                            sin(bearing(j)), range(j)*cos(bearing(j))];
                   Jac(2*meas_no -1:2*meas_no,2*i-1:2*i) = [-cos(bearing(i)), range(i)*sin(bearing(i));...
                                                            -sin(bearing(i)), -range(i)*cos(bearing(i))];
                   meas_no = meas_no + 1;                                   
               end

            end

        end
    end
    R = [R_range; R_bearing;];
     
    covariance = Jac*diag(repmat(R,17,1))*Jac.';
    covariance = covariance + eps.*(eye(length(covariance)));
    

end

