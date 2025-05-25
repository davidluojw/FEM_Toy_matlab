function [k_ele, f_ele] = elasticity_elemTri6_2d(model, ee)

    % =========================================================================
    % generate the quadrature rule
    [xi, eta, zeta, weight] = gauss_Tri3_2d(model.ngp);
    % =========================================================================
    k_ele = zeros(model.nee, model.nee);   % dimension of element stiffness is 6 x 6
    f_ele = zeros(model.nee, 1);           % element force nodes

    nodes_ele = zeros(model.nen, model.ndof);  % the nodes in the ee-th element, (x1,y1;x2,y2,x3,y3)
    
    for aa = 1 : model.nen
        for ii = 1: model.ndof 
            nodes_ele(aa, ii) = model.nodes( model.IEN(aa,ee), ii);
        end  
    end

    % loop over quadrature points   
    for ll = 1 : model.nint

        x_l = 0.0; y_l = 0.0;           % coordinate in terms of xi(ll)
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;

        % shape function matrix
        N1 = T6_ShapeBasisN_2d(1, xi(ll), eta(ll), zeta(ll));
        N2 = T6_ShapeBasisN_2d(2, xi(ll), eta(ll), zeta(ll));
        N3 = T6_ShapeBasisN_2d(3, xi(ll), eta(ll), zeta(ll));
        N4 = T6_ShapeBasisN_2d(4, xi(ll), eta(ll), zeta(ll));
        N5 = T6_ShapeBasisN_2d(5, xi(ll), eta(ll), zeta(ll));
        N6 = T6_ShapeBasisN_2d(6, xi(ll), eta(ll), zeta(ll));

        x_l = [N1, N2, N3, N4, N5, N6] * nodes_ele(:, 1);
        y_l = [N1, N2, N3, N4, N5, N6] * nodes_ele(:, 2);

        N = [N1, 0.0, N2, 0.0, N3, 0.0, N4, 0.0, N5, 0.0, N6, 0.0;
             0.0, N1, 0.0, N2, 0.0, N3, 0.0, N4, 0.0, N5, 0.0, N6];
        
        % Grad of shape function matrix, BB  = (N1_xi,  N2_xi,  N3_xi, N4_xi,  N5_xi,  N6_xi; 
        %                                       N1_eta, N2_eta, N3_eta, N4_eta, N5_eta, N6_eta)
        [N1_xi, N1_eta] = T6_ShapeBasisN_Grad_2d(1, xi(ll), eta(ll));
        [N2_xi, N2_eta] = T6_ShapeBasisN_Grad_2d(2, xi(ll), eta(ll));
        [N3_xi, N3_eta] = T6_ShapeBasisN_Grad_2d(3, xi(ll), eta(ll));
        [N4_xi, N4_eta] = T6_ShapeBasisN_Grad_2d(4, xi(ll), eta(ll));
        [N5_xi, N5_eta] = T6_ShapeBasisN_Grad_2d(5, xi(ll), eta(ll));
        [N6_xi, N6_eta] = T6_ShapeBasisN_Grad_2d(6, xi(ll), eta(ll));

        GradN = [N1_xi,  N2_xi,  N3_xi,  N4_xi,  N5_xi,  N6_xi;
                 N1_eta, N2_eta, N3_eta, N4_eta, N5_eta, N6_eta];
        
        % Jacobian matrix
        J = GradN * nodes_ele;
        detJ = det(J);

        % Grad of shape function matrix, BB  = (N1_x, N2_x, N3_x, N4_x, N5_x, N6_x; 
        %                                       N1_y, N2_y, N3_y, N4_y, N5_y, N6_y)
        % J * BB = GradN, x_xi * N1_x + y_xi * N1_y = N1_xi
        %                 x_xi * N2_x + y_xi * N2_y = N2_xi
        %                 x_xi * N3_x + y_xi * N3_y = N3_xi
        %                 x_xi * N4_x + y_xi * N4_y = N4_xi
        %                 x_xi * N5_x + y_xi * N5_y = N5_xi
        %                 x_xi * N6_x + y_xi * N6_y = N6_xi
        %                 x_eta * N1_x + y_eta * N1_y = N1_eta
        %                 x_eta * N2_x + y_eta * N2_y = N2_eta
        %                 x_eta * N3_x + y_eta * N3_y = N3_eta
        %                 x_eta * N4_x + y_eta * N4_y = N4_eta
        %                 x_eta * N5_x + y_eta * N5_y = N5_eta
        %                 x_eta * N6_x + y_eta * N6_y = N6_eta
        BB = J \ GradN;
        B1_x = BB(1,1); B2_x = BB(1,2); B3_x = BB(1,3); 
        B4_x = BB(1,4); B5_x = BB(1,5); B6_x = BB(1,6);
        B1_y = BB(2,1); B2_y = BB(2,2); B3_y = BB(2,3);
        B4_y = BB(2,4); B5_y = BB(2,5); B6_y = BB(2,6);

        Bmat = [B1_x, 0.0,  B2_x, 0.0,  B3_x, 0.0, B4_x, 0.0,  B5_x, 0.0,  B6_x, 0.0;
                0.0,  B1_y, 0.0,  B2_y, 0.0,  B3_y, 0.0,  B4_y, 0.0,  B5_y, 0.0,  B6_y;
                B1_y, B1_x, B2_y, B2_x, B3_y, B3_x, B4_y, B4_x, B5_y, B5_x, B6_y, B6_x];
        
        % element stiffness maatrix
        k_ele = k_ele + Bmat' * model.D * Bmat * detJ * weight(ll);

        % element nodal forces, body forces
        be = [model.exact_bf_x(x_l, y_l); model.exact_bf_y(x_l, y_l)];
        f_ele = f_ele + N' * be * detJ * weight(ll);

    end % end of quadrature loop

        
end
