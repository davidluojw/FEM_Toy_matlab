function model = get_nodal_stress(model, ee)
    % NODAL_STRESS element nodal stress and cumulate to global stress

    for ee = 1:model.nel
        % element nodal displacement
        de = model.d(model.LM(:, ee)); 

        % element node coordinates
        nodes_ele_idx = model.IEN(:, ee);
        nodes_ele = model.nodes(nodes_ele_idx, :);
        
        % xi-eta-zeta coordinate evaluation, nodal points
        xi = [1, 0, 0];
        eta = [0, 1, 0];
        strain = zeros(3, model.nen);
        stress = zeros(3, model.nen);
        
        % loop over nodal points
        for ii = 1 : model.nen
            
            % Grad of shape function matrix, BB  = (N1_xi,  N2_xi,  N3_xi; 
            %                                       N1_eta, N2_eta, N3_eta)
            [N1_xi, N1_eta] = Area_ShapeBasisN_Grad_2d(model.degree, 1, xi(ii), eta(ii));
            [N2_xi, N2_eta] = Area_ShapeBasisN_Grad_2d(model.degree, 2, xi(ii), eta(ii));
            [N3_xi, N3_eta] = Area_ShapeBasisN_Grad_2d(model.degree, 3, xi(ii), eta(ii));

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