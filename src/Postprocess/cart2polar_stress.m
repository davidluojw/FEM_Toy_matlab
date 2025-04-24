function  cart2polar_stress(model)
    %  transform stress from cartesian coordiante to polar coordinate

    fprintf('\ts_r\t\ts_theta\n');

    for ii = 1:model.nnp

        XX = model.nodes(ii, 1);
        YY = model.nodes(ii, 2);

        theta = atan2(YY, XX);

        cos2theta = cos(2*theta);
        sin2theta = sin(2*theta);

        sigma_xx = model.stress_nodal(1, ii);
        sigma_yy = model.stress_nodal(2, ii);
        sigma_xy = model.stress_nodal(3, ii);


        sigma_r = 0.5*(sigma_xx + sigma_yy) + 0.5*(sigma_xx - sigma_yy).*cos2theta + ...
                  sigma_xy.*sin2theta;
    
        sigma_theta = 0.5*(sigma_xx + sigma_yy) - 0.5*(sigma_xx - sigma_yy).*cos2theta - ...
                    sigma_xy.*sin2theta;

        fprintf('\t%.6e\t%.6e\n', sigma_r, sigma_theta);
        
    end



end
