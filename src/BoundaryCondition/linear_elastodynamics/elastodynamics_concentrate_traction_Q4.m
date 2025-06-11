function model = elastodynamics_concentrate_traction_Q4(model, tt)

    % concentrate force
    for ii = 1:model.neq
        model.f(ii) = model.f(ii) + model.P(ii);
    end

    for ii = 1:size(model.nbc_nodes,2)
        BCedge = model.nbc_nodes(:, ii);

        if BCedge(1) == 0
            break;
        end

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

        [qpg1D, weight1D] = gauss(model.ngp, -1, 1);

        for ll = 1:model.ngp
            % boudary edge simplied shape function
            N1 = Poly_ShapeBasisN_1d(1, 1, qpg1D(ll), 0);
            N2 = Poly_ShapeBasisN_1d(1, 2, qpg1D(ll), 0);

            Nt = [N1, 0, N2, 0;...
                  0, N1, 0, N2];

            % the position of the quadrature points
            x_l_h = [N1, N2] * x_ele;
            y_l_h = [N1, N2] * y_ele;

            % traction force
            cur_t = model.time_intervals(tt); 
            t_l = [model.exact_t_x(x_l_h, y_l_h, cur_t); model.exact_t_y(x_l_h, y_l_h, cur_t)];

            model.f(dof_ids) = model.f(dof_ids) + weight1D(ll) * detJ * Nt' * t_l;
        end
    end

end
