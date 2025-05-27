function  model =  error_analysis_Q8(model)

    % postprocess the solution by calculating the error measured in L2 norm
    nqp = 10; % we need more points 
    [xi, eta, weight] = gauss_2d(model.ngp, model.ngp);

    % Extract the displacement component
    model.ux_h = zeros(model.neq / 2, 1);
    model.uy_h = zeros(model.neq / 2, 1);
    for ii = 1 : model.neq
        if mod(ii, 2) == 1
            model.ux_h(floor(ii / 2) + 1) = model.d(ii);
        else
            model.uy_h(floor(ii / 2)) = model.d(ii);
        end
    end
    
    model.errorL2_x = 0.0; model.errorL2_y = 0.0; bottomL2_x = 0.0; bottomL2_y = 0.0;
    model.errorH1_x = 0.0; model.errorH1_y = 0.0; bottomH1_x = 0.0; bottomH1_y = 0.0;
    for ee = 1 : model.nel
        ele_ids = model.IEN(1:model.nen, ee);  % the number of nodes within the ee element
        x_ele = model.nodes(ele_ids, 1); 
        y_ele = model.nodes(ele_ids, 2);

        ux_ele = model.ux_h( ele_ids );
        uy_ele = model.uy_h( ele_ids );

        for ll = 1:model.nint
            x_l = 0.0; y_l = 0.0;    % the coordinates at the quadrature points
            ux_l = 0.0; uy_l = 0.0;  % the displacements at the quadrature points
            ux_l_xi = 0.0; ux_l_eta = 0.0;   % the derivative of the displacement at the quadrature points
            uy_l_xi = 0.0; uy_l_eta = 0.0;
            dx_dxi = 0.0; dy_dxi = 0.0;    % the derivative of the coordinates at the quadrature points
            dx_deta = 0.0; dy_deta = 0.0;

            for aa = 1:model.nen

                x_l = x_l + x_ele(aa) * Q8_ShapeBasisN_2d(aa, xi(ll), eta(ll));
                y_l = y_l + y_ele(aa) * Q8_ShapeBasisN_2d(aa, xi(ll), eta(ll));
                ux_l = ux_l + ux_ele(aa) * Q8_ShapeBasisN_2d(aa, xi(ll), eta(ll));
                uy_l = uy_l + uy_ele(aa) * Q8_ShapeBasisN_2d(aa, xi(ll), eta(ll));

                [Na_xi, Na_eta] = Q8_ShapeBasisN_Grad_2d(aa, xi(ll), eta(ll));
                ux_l_xi  = ux_l_xi  + ux_ele(aa) * Na_xi;
                uy_l_xi  = uy_l_xi  + uy_ele(aa) * Na_xi;
                ux_l_eta = ux_l_eta + ux_ele(aa) * Na_eta;
                uy_l_eta = uy_l_eta + uy_ele(aa) * Na_eta;
                dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                dx_deta = dx_deta + x_ele(aa) * Na_eta;
                dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                dy_deta = dy_deta + y_ele(aa) * Na_eta;

            end
            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    
            ux_l_x = (ux_l_xi * dy_deta - ux_l_eta * dy_dxi) / detJ;
            ux_l_y = (ux_l_xi * (-dx_deta) + ux_l_eta * dx_dxi) / detJ;

            uy_l_x = (uy_l_xi * dy_deta - uy_l_eta * dy_dxi) / detJ;
            uy_l_y = (uy_l_xi * (-dx_deta) + uy_l_eta * dx_dxi) / detJ;
            
            model.errorL2_x = model.errorL2_x + weight(ll) * detJ * (ux_l - model.exact_ux(x_l, y_l))^2;
            model.errorL2_y = model.errorL2_y + weight(ll) * detJ * (uy_l - model.exact_uy(x_l, y_l))^2;

            model.errorH1_x = model.errorH1_x + weight(ll) * detJ * ( ( ux_l_x- model.exact_ux_x(x_l,y_l))^2 + ( ux_l_y - model.exact_ux_y(x_l,y_l))^2 );
            model.errorH1_y = model.errorH1_y + weight(ll) * detJ * ( ( uy_l_x- model.exact_uy_x(x_l,y_l))^2 + ( uy_l_y - model.exact_uy_y(x_l,y_l))^2 );

            bottomL2_x = bottomL2_x + weight(ll) * detJ * model.exact_ux(x_l, y_l)^2;
            bottomL2_y = bottomL2_y + weight(ll) * detJ * model.exact_uy(x_l, y_l)^2;

            bottomH1_x = bottomH1_x + weight(ll) * detJ * ( model.exact_ux_x(x_l,y_l)^2 + model.exact_ux_y(x_l,y_l)^2 );
            bottomH1_y = bottomH1_y + weight(ll) * detJ * ( model.exact_uy_x(x_l,y_l)^2 + model.exact_uy_y(x_l,y_l)^2 );
        end
    end
    
    model.errorL2_x = sqrt(model.errorL2_x) / sqrt(bottomL2_x);
    model.errorL2_y = sqrt(model.errorL2_y) / sqrt(bottomL2_y);

    model.errorH1_x = sqrt(model.errorH1_x) / sqrt(bottomH1_x);
    model.errorH1_y = sqrt(model.errorH1_y) / sqrt(bottomH1_y);
    
    model.errorL2 = model.errorL2_x + model.errorL2_y;
    model.errorH1 = model.errorH1_x + model.errorH1_y;


end

