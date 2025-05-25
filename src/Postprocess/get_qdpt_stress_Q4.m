function model = get_qdpt_stress_Q4(model)
    % GET_STRESS: QUADRATURE POINT STRESS
    
    for ee =1:model.nel
        % element nodal displacement
        de = model.d(model.LM(:, ee)); 
        
        % element node coordinates
        nodes_ele = model.nodes(model.IEN(:, ee), :);
        
        % quadrature points
        [xi, eta, weight] = gauss_2d(model.ngp, model.ngp);
        
        % initialization
        nodes_l = zeros(model.nint, 2);
        
        % loop over quadrature points   
        for ll = 1 : model.nint

            % shape function matrix
            N1 = Quad_ShapeBasisN_2d(model.deg_xi, model.deg_eta, 1, xi(ll), eta(ll));
            N2 = Quad_ShapeBasisN_2d(model.deg_xi, model.deg_eta, 2, xi(ll), eta(ll));
            N3 = Quad_ShapeBasisN_2d(model.deg_xi, model.deg_eta, 3, xi(ll), eta(ll));
            N4 = Quad_ShapeBasisN_2d(model.deg_xi, model.deg_eta, 4, xi(ll), eta(ll));

            N = [N1, 0.0, N2, 0.0, N3, 0.0, N4, 0.0;
                0.0, N1, 0.0, N2, 0.0, N3, 0.0, N4];
            
            % Grad of shape function matrix, BB  = (N1_xi,  N2_xi,  N3_xi,  N4_xi; 
            %                                       N1_eta, N2_eta, N3_eta, N4_eta)
            [N1_xi, N1_eta] = Quad_ShapeBasisN_Grad_2d(model.deg_xi, model.deg_eta, 1, xi(ll), eta(ll));
            [N2_xi, N2_eta] = Quad_ShapeBasisN_Grad_2d(model.deg_xi, model.deg_eta, 2, xi(ll), eta(ll));
            [N3_xi, N3_eta] = Quad_ShapeBasisN_Grad_2d(model.deg_xi, model.deg_eta, 3, xi(ll), eta(ll));
            [N4_xi, N4_eta] = Quad_ShapeBasisN_Grad_2d(model.deg_xi, model.deg_eta, 4, xi(ll), eta(ll));

            GradN = [N1_xi,  N2_xi,  N3_xi, N4_xi;
                    N1_eta, N2_eta, N3_eta, N4_eta];
            
            % Jacobian matrix
            J = GradN * nodes_ele;
            detJ = det(J);

            % Grad of shape function matrix, BB  = (N1_x, N2_x, N3_x, N4_x; 
            %                                       N1_y, N2_y, N3_y, N4_y)
            % J * BB = GradN, x_xi * N1_x + y_xi * N1_y = N1_xi
            %                 x_xi * N2_x + y_xi * N2_y = N2_xi
            %                 x_xi * N3_x + y_xi * N3_y = N3_xi
            %                 x_xi * N4_x + y_xi * N4_y = N4_xi
            %                 x_eta * N1_x + y_eta * N1_y = N1_eta
            %                 x_eta * N2_x + y_eta * N2_y = N2_eta
            %                 x_eta * N3_x + y_eta * N3_y = N3_eta
            %                 x_eta * N4_x + y_eta * N4_y = N4_eta
            BB = J \ GradN;
            B1_x = BB(1,1); B2_x = BB(1,2); B3_x = BB(1,3); B4_x = BB(1,4);
            B1_y = BB(2,1); B2_y = BB(2,2); B3_y = BB(2,3); B4_y = BB(2,4);

            Bmat = [B1_x, 0.0,  B2_x, 0.0,  B3_x, 0.0,  B4_x, 0.0;
                    0.0,  B1_y, 0.0,  B2_y, 0.0,  B3_y, 0.0,  B4_y;
                    B1_y, B1_x, B2_y, B2_x, B3_y, B3_x, B4_y, B4_x];
            
            % Na, quadrature point coordinates
            Na = [N1, N2, N3, N4];
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