function model = concentrate_traction_Q8(model)

    % concentrate force
    for ii = 1:model.neq
        model.f(ii) = model.f(ii) + model.P(ii);
    end

    [xi, weight] = gauss(model.ngp, -1, 1);

    % natural bc: traction force (edge)
    for ii = 1:model.nbe
        ft = zeros(6, 1);
        node1 = model.n_bc(1, ii);
        node2 = model.n_bc(2, ii);
        node3 = model.n_bc(3, ii);

        n_bce = model.n_bc(4:end, ii);

        x1 = model.nodes(node1, 1);
        y1 = model.nodes(node1, 2);
        x2 = model.nodes(node2, 1);
        y2 = model.nodes(node2, 2);
        x3 = model.nodes(node3, 1);
        y3 = model.nodes(node3, 2);


        leng = sqrt((x3-x1)^2 + (y3-y1)^2);
        J = leng/ 2.0;

        for ll = 1:model.ngp
            % shape function matrix
            N = 0.5*[xi(ll)*(xi(ll) - 1), 0, 2*(1-xi(ll)*xi(ll)), 0, xi(ll)*(xi(ll) + 1), 0;   
                     0, xi(ll)*(xi(ll) - 1), 0, 2*(1-xi(ll)*xi(ll)), 0, xi(ll)*(xi(ll) + 1)];
            NQ = 0.5*[1-xi(ll), 0, 1+xi(ll), 0;   
                      0, 1-xi(ll), 0, 1+xi(ll)];
            
            traction  = NQ * n_bce;
            ft = ft + weight(ll) * J * N' * traction;
        end

        model.f(model.ID(1, node1)) = model.f(model.ID(1, node1)) + ft(1);
        model.f(model.ID(2, node1)) = model.f(model.ID(2, node1)) + ft(2);
        model.f(model.ID(1, node2)) = model.f(model.ID(1, node2)) + ft(3);
        model.f(model.ID(2, node2)) = model.f(model.ID(2, node2)) + ft(4);
        model.f(model.ID(1, node3)) = model.f(model.ID(1, node3)) + ft(5);
        model.f(model.ID(2, node3)) = model.f(model.ID(2, node3)) + ft(6);

    end
    

end
