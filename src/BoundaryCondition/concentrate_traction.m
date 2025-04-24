function model = concentrate_traction(model)

    % concentrate force
    for ii = 1:model.neq
        model.f(ii) = model.f(ii) + model.P(ii);
    end

    [xi, weight] = gauss(model.nint, -1, 1);

    % natural bc: traction force (edge)
    for ii = 1:model.nbe
        ft = zeros(4, 1);
        node1 = model.n_bc(1, ii);
        node2 = model.n_bc(2, ii);

        n_bce = model.n_bc(3:end, ii);

        x1 = model.nodes(node1, 1);
        y1 = model.nodes(node1, 2);
        x2 = model.nodes(node2, 1);
        y2 = model.nodes(node2, 2);

        leng = sqrt((x2-x1)^2 + (y2-y1)^2);
        J = leng/ 2.0;

        for ll = 1:model.nint
            N = 0.5*[1-xi(ll), 0, 1+xi(ll), 0;
                     0, 1-xi(ll), 0, 1+xi(ll)];
            
            traction  = N * n_bce;
            ft = ft + weight(ll) * J * N' * traction;
        end

        model.f(model.ID(1, node1)) = model.f(model.ID(1, node1)) + ft(1);
        model.f(model.ID(2, node1)) = model.f(model.ID(2, node1)) + ft(2);
        model.f(model.ID(1, node2)) = model.f(model.ID(1, node2)) + ft(3);
        model.f(model.ID(2, node2)) = model.f(model.ID(2, node2)) + ft(4);

    end
    

end
