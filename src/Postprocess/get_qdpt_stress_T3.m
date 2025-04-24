function model = get_qdpt_stress_T3(model)
    % GET_STRESS: QUADRATURE POINT STRESS
    
    for ee =1:model.nel
        % element nodal displacement
        de = model.d(model.LM(:, ee)); 
        
        % element node coordinates
        nodes_ele = model.nodes(model.IEN(:, ee), :);
        
        % quadrature points
        [xi, eta, zeta, weight] = gauss_Tri3_2d(model.ngp);
        
        % initialization
        nodes_l = zeros(model.nint, 2);
        
        % loop over quadrature points   
        for ll = 1 : model.nint

            % shape function matrix
            N1 = Area_ShapeBasisN_2d(model.degree, 1, xi(ll), eta(ll), zeta(ll));
            N2 = Area_ShapeBasisN_2d(model.degree, 2, xi(ll), eta(ll), zeta(ll));
            N3 = Area_ShapeBasisN_2d(model.degree, 3, xi(ll), eta(ll), zeta(ll));

            N = [N1, 0.0, N2, 0.0, N3, 0.0;
                0.0, N1, 0.0, N2, 0.0, N3];
            
            % Grad of shape function matrix, BB  = (N1_xi,  N2_xi,  N3_xi; 
            %                                       N1_eta, N2_eta, N3_eta)
            [N1_xi, N1_eta] = Area_ShapeBasisN_Grad_2d(model.degree, 1, xi(ll), eta(ll));
            [N2_xi, N2_eta] = Area_ShapeBasisN_Grad_2d(model.degree, 2, xi(ll), eta(ll));
            [N3_xi, N3_eta] = Area_ShapeBasisN_Grad_2d(model.degree, 3, xi(ll), eta(ll));

            GradN = [N1_xi,  N2_xi,  N3_xi;
                    N1_eta, N2_eta, N3_eta];
            
            % Jacobian matrix
            J = GradN * nodes_ele;
            detJ = det(J);

            % Grad of shape function matrix, BB  = (N1_x, N2_x, N3_x; 
            %                                       N1_y, N2_y, N3_y)
            % J * BB = GradN, x_xi * N1_x + y_xi * N1_y = N1_xi
            %                 x_xi * N2_x + y_xi * N2_y = N2_xi
            %                 x_xi * N3_x + y_xi * N3_y = N3_xi
            %                 x_eta * N1_x + y_eta * N1_y = N1_eta
            %                 x_eta * N2_x + y_eta * N2_y = N2_eta
            %                 x_eta * N3_x + y_eta * N3_y = N3_eta
            BB = J \ GradN;
            B1_x = BB(1,1); B2_x = BB(1,2); B3_x = BB(1,3); 
            B1_y = BB(2,1); B2_y = BB(2,2); B3_y = BB(2,3);

            Bmat = [B1_x, 0.0,  B2_x, 0.0,  B3_x, 0.0;
                    0.0,  B1_y, 0.0,  B2_y, 0.0,  B3_y;
                    B1_y, B1_x, B2_y, B2_x, B3_y, B3_x];
            
            % Na, quadrature point coordinates
            Na = [N1, N2, N3];
            nodes_l(ll, :) = Na * nodes_ele; 

            % strain
            model.strain_qdpt(:, ll, ee) = Bmat * de;

            % stress
            model.stress_qdpt(:, ll, ee) = model.D * model.strain_qdpt(:, ll, ee);
        end
        
        % print results
        fprintf('\tx-coord\t\ty-coord\t\ts_xx\t\ts_yy\t\ts_xy\n');
        for ll = 1:model.nint
            fprintf('\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n', nodes_l(ll, 1), nodes_l(ll, 2), ...
                                                        model.stress_qdpt(1, ll, ee), model.stress_qdpt(2, ll, ee), model.stress_qdpt(3, ll, ee));
        end
    end
end