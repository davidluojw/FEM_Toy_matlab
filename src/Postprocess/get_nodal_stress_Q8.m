function model = get_nodal_stress_Q8(model, ee)
    % NODAL_STRESS element nodal stress and cumulate to global stress

    for ee = 1:model.nel
        % element nodal displacement
        de = model.d(model.LM(:, ee)); 

        % element node coordinates
        nodes_ele_idx = model.IEN(:, ee);
        nodes_ele = model.nodes(nodes_ele_idx, :);
        
        % xi-eta coordinate evaluation, nodal points
        xi =  [-1,  1, 1, -1, 0,  1, 0, -1];
        eta = [-1, -1, 1,  1, -1, 0, 1, 0];
        strain = zeros(3, model.nen);
        stress = zeros(3, model.nen);
        
        % loop over nodal points
        for ii = 1 : model.nen
            
            % Grad of shape function matrix, BB  = (N1_xi,  N2_xi,  N3_xi,  N4_xi, N5_xi, N6_xi, N7_xi, N8_xi; 
            %                                       N1_eta, N2_eta, N3_eta, N4_eta, N5_eta, N6_eta, N7_eta, N8_eta)
            [N1_xi, N1_eta] = Q8_ShapeBasisN_Grad_2d(1, xi(ii), eta(ii));
            [N2_xi, N2_eta] = Q8_ShapeBasisN_Grad_2d(2, xi(ii), eta(ii));
            [N3_xi, N3_eta] = Q8_ShapeBasisN_Grad_2d(3, xi(ii), eta(ii));
            [N4_xi, N4_eta] = Q8_ShapeBasisN_Grad_2d(4, xi(ii), eta(ii));
            [N5_xi, N5_eta] = Q8_ShapeBasisN_Grad_2d(5, xi(ii), eta(ii));
            [N6_xi, N6_eta] = Q8_ShapeBasisN_Grad_2d(6, xi(ii), eta(ii));
            [N7_xi, N7_eta] = Q8_ShapeBasisN_Grad_2d(7, xi(ii), eta(ii));
            [N8_xi, N8_eta] = Q8_ShapeBasisN_Grad_2d(8, xi(ii), eta(ii));

            GradN = [N1_xi,  N2_xi,  N3_xi,  N4_xi,  N5_xi,  N6_xi,  N7_xi,  N8_xi;
                     N1_eta, N2_eta, N3_eta, N4_eta, N5_eta, N6_eta, N7_eta, N8_eta];
            
            % Jacobian matrix
            J = GradN * nodes_ele;
            detJ = det(J);

            % Grad of shape function matrix, BB  = (N1_x, N2_x, N3_x, N4_x, N5_x, N6_x, N7_x, N8_x; 
            %                                       N1_y, N2_y, N3_y, N4_y, N6_y, N7_y, N8_y, N9_y)
            % J * BB = GradN, x_xi * N1_x + y_xi * N1_y = N1_xi
            %                 x_xi * N2_x + y_xi * N2_y = N2_xi
            %                 x_xi * N3_x + y_xi * N3_y = N3_xi
            %                 x_xi * N4_x + y_xi * N4_y = N4_xi
            %                 x_xi * N5_x + y_xi * N5_y = N5_xi
            %                 x_xi * N6_x + y_xi * N6_y = N6_xi
            %                 x_xi * N7_x + y_xi * N7_y = N7_xi
            %                 x_xi * N8_x + y_xi * N8_y = N8_xi
            %                 x_eta * N1_x + y_eta * N1_y = N1_eta
            %                 x_eta * N2_x + y_eta * N2_y = N2_eta
            %                 x_eta * N3_x + y_eta * N3_y = N3_eta
            %                 x_eta * N4_x + y_eta * N4_y = N4_eta
            %                 x_eta * N5_x + y_eta * N5_y = N5_eta
            %                 x_eta * N6_x + y_eta * N6_y = N6_eta
            %                 x_eta * N7_x + y_eta * N7_y = N7_eta
            %                 x_eta * N8_x + y_eta * N8_y = N8_eta
            BB = J \ GradN;
            B1_x = BB(1,1); B2_x = BB(1,2); B3_x = BB(1,3); B4_x = BB(1,4);
            B5_x = BB(1,5); B6_x = BB(1,6); B7_x = BB(1,7); B8_x = BB(1,8);
            B1_y = BB(2,1); B2_y = BB(2,2); B3_y = BB(2,3); B4_y = BB(2,4);
            B5_y = BB(2,5); B6_y = BB(2,6); B7_y = BB(2,7); B8_y = BB(2,8);

            Bmat = [B1_x, 0.0,  B2_x, 0.0,  B3_x, 0.0,  B4_x, 0.0,  B5_x, 0.0,  B6_x, 0.0,  B7_x, 0.0,  B8_x, 0.0;
                    0.0,  B1_y, 0.0,  B2_y, 0.0,  B3_y, 0.0,  B4_y, 0.0,  B5_y, 0.0,  B6_y, 0.0,  B7_y, 0.0,  B8_y;
                    B1_y, B1_x, B2_y, B2_x, B3_y, B3_x, B4_y, B4_x, B5_y, B5_x, B6_y, B6_x, B7_y, B7_x, B8_y, B8_x];
            

            % strain
            strain(:, ii) = Bmat * de;

            % stress
            stress(:, ii) = model.D * strain(:, ii);
        end

        
        % accumulate the counter of stress and value
        model.counter_adjpt(nodes_ele_idx) = model.counter_adjpt(nodes_ele_idx) + ones(1, model.nen);
        model.strain_nodal(:, nodes_ele_idx) = model.strain_nodal(:, nodes_ele_idx) + strain;
        model.stress_nodal(:, nodes_ele_idx) = model.stress_nodal(:, nodes_ele_idx) + stress;
    end

    % mean value of the adjent element stress
    model.strain_nodal = model.strain_nodal ./ model.counter_adjpt;
    model.stress_nodal = model.stress_nodal ./ model.counter_adjpt;

    % print results
    fprintf('\tx-coord\t\ty-coord\t\ts_xx\t\ts_yy\t\ts_xy\n');
    for ii = 1:model.nnp
        fprintf('\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n', model.nodes(ii, 1), model.nodes(ii, 2), ...
                model.stress_nodal(1, ii), model.stress_nodal(2, ii), model.stress_nodal(3, ii));
    end
end