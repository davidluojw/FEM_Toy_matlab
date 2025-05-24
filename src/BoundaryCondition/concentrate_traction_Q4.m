function model = concentrate_traction_Q4(model)

    % concentrate force
    for ii = 1:model.neq
        model.f(ii) = model.f(ii) + model.P(ii);
    end

    % natural bc: traction force (edge)
    % for ii = 1:model.nbe
    %     ft = zeros(4, 1);
    %     node1 = model.n_bc(1, ii);
    %     node2 = model.n_bc(2, ii);

    %     n_bce = model.n_bc(3:end, ii);

    %     x1 = model.nodes(node1, 1);
    %     y1 = model.nodes(node1, 2);
    %     x2 = model.nodes(node2, 1);
    %     y2 = model.nodes(node2, 2);

    %     leng = sqrt((x2-x1)^2 + (y2-y1)^2);
    %     J = leng/ 2.0;

    %     for ll = 1:model.ngp
    %         % shape function matrix
    %         N = 0.5*[1-xi1D(ll), 0, 1+xi1D(ll), 0;   
    %                  0, 1-xi1D(ll), 0, 1+xi1D(ll)];
            
    %         traction  = N * n_bce;
    %         ft = ft + weight1D(ll) * J * N' * traction;
    %     end

    %     model.f(model.ID(1, node1)) = model.f(model.ID(1, node1)) + ft(1);
    %     model.f(model.ID(2, node1)) = model.f(model.ID(2, node1)) + ft(2);
    %     model.f(model.ID(1, node2)) = model.f(model.ID(1, node2)) + ft(3);
    %     model.f(model.ID(2, node2)) = model.f(model.ID(2, node2)) + ft(4);

    % end


    for ii = 1:size(model.nbc_nodes,2)
        BCedge = model.nbc_nodes(:, ii);

        node1 = BCedge(1);
        node2 = BCedge(2);

        dof_ids = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2]';

        x1 = model.nodes(node1, 1);
        y1 = model.nodes(node1, 2);
        x2 = model.nodes(node2, 1);
        y2 = model.nodes(node2, 2);
        len = sqrt((x2-x1)^2 + (y2-y1)^2);
        detJ = len / 2;

        x_ele = [x1; x2];
        y_ele = [y1; y2];

        [qpg1D, weight1D] = gauss(3, -1, 1);

        for ll = 1:3
            N1 = Poly_ShapeBasisN_1d(1, 1, qpg1D(ll), 0);
            N2 = Poly_ShapeBasisN_1d(1, 2, qpg1D(ll), 0);

            % shape function matrix
            % N1 = Quad_ShpaeBasisN_2d(model.deg_xi, model.deg_eta, 4, qpg1D(ll), 1);
            % N2 = Quad_ShpaeBasisN_2d(model.deg_xi, model.deg_eta, 3, qpg1D(ll), 1);
            % N3 = Quad_ShpaeBasisN_2d(model.deg_xi, model.deg_eta, 3, xi1D(ll), -1);  % N3 = 0
            % N4 = Quad_ShpaeBasisN_2d(model.deg_xi, model.deg_eta, 4, xi1D(ll), -1);  % N4 = 0

            Nt = [N1, 0, N2, 0;...
                  0, N1, 0, N2];

            x_l_h = [N1, N2] * x_ele;
            y_l_h = [N1, N2] * y_ele;


            t_l = [model.exact_t_x(x_l_h, y_l_h); model.exact_t_y(x_l_h, y_l_h)];

            model.f(dof_ids) = model.f(dof_ids) + weight1D(ll) * detJ * Nt' * t_l;
        end
    end

end
