function model = concentrate_traction_Q8(model)

    % concentrate force
    for ii = 1:model.neq
        model.f(ii) = model.f(ii) + model.P(ii);
    end


    for ii = 1:size(model.nbc_nodes,2)
        BCedge = model.nbc_nodes(:, ii);

        node1 = BCedge(1);
        node2 = BCedge(2);
        node3 = BCedge(3);

        dof_ids = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2, 2*node3 - 1, 2*node3]';

        x1 = model.nodes(node1, 1);
        y1 = model.nodes(node1, 2);
        x2 = model.nodes(node2, 1);
        y2 = model.nodes(node2, 2);
        x3 = model.nodes(node3, 1);
        y3 = model.nodes(node3, 2);
        len = sqrt((x3-x1)^2 + (y3-y1)^2);
        detJ = len / 2;

        x_ele = [x1; x2; x3];
        y_ele = [y1; y2; y3];

        [qpg1D, weight1D] = gauss(model.ngp, -1, 1);

        for ll = 1:model.ngp
            % boudary edge simplied shape function
            N1 = Poly_ShapeBasisN_1d(2, 1, qpg1D(ll), 0);
            N2 = Poly_ShapeBasisN_1d(2, 2, qpg1D(ll), 0);
            N3 = Poly_ShapeBasisN_1d(2, 3, qpg1D(ll), 0);

            Nt = [N1, 0, N2, 0, N3, 0;...
                  0, N1, 0, N2, 0, N3];

            % the position of the quadrature points
            x_l_h = [N1, N2, N3] * x_ele;
            y_l_h = [N1, N2, N3] * y_ele;

            % traction force
            t_l = [model.exact_t_x(x_l_h, y_l_h); model.exact_t_y(x_l_h, y_l_h)];

            model.f(dof_ids) = model.f(dof_ids) + weight1D(ll) * detJ * Nt' * t_l;
        end
    end

    

end
