function [k_ele, f_ele] = elasticity_elemQuad4_2d(model, ee)

    % =========================================================================
    % generate the quadrature rule
    [xi, eta, weight] = gauss_2d(model.ngp, model.ngp);
    [xi1D, weight1D] = gauss(model.ngp, -1, 1);
    % =========================================================================
    k_ele = zeros(model.nee, model.nee);   % dimension of element stiffness is 8 x 8
    f_ele = zeros(model.nee, 1);           % element force nodes

    nodes_ele = zeros(model.nen, model.ndof);  % the nodes in the ee-th element, (x1,y1;x2,y2,x3,y3,x4,y4)
    
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
        N1 = Quad_ShpaeBasisN_2d(model.deg_xi, model.deg_eta, 1, xi(ll), eta(ll));
        N2 = Quad_ShpaeBasisN_2d(model.deg_xi, model.deg_eta, 2, xi(ll), eta(ll));
        N3 = Quad_ShpaeBasisN_2d(model.deg_xi, model.deg_eta, 3, xi(ll), eta(ll));
        N4 = Quad_ShpaeBasisN_2d(model.deg_xi, model.deg_eta, 4, xi(ll), eta(ll));

        N = [N1, 0.0, N2, 0.0, N3, 0.0, N4, 0.0;
             0.0, N1, 0.0, N2, 0.0, N3, 0.0, N4];
        
        % Grad of shape function matrix, BB  = (N1_xi,  N2_xi,  N3_xi,  N4_xi; 
        %                                       N1_eta, N2_eta, N3_eta, N4_eta)
        [N1_xi, N1_eta] = Quad_ShapeBasisN_Grad_2d(model.deg_xi, model.deg_eta, 1, xi(ll), eta(ll));
        [N2_xi, N2_eta] = Quad_ShapeBasisN_Grad_2d(model.deg_xi, model.deg_eta, 2, xi(ll), eta(ll));
        [N3_xi, N3_eta] = Quad_ShapeBasisN_Grad_2d(model.deg_xi, model.deg_eta, 3, xi(ll), eta(ll));
        [N4_xi, N4_eta] = Quad_ShapeBasisN_Grad_2d(model.deg_xi, model.deg_eta, 4, xi(ll), eta(ll));

        GradN = [N1_xi,  N2_xi,  N3_xi,  N4_xi;
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
        
        % element stiffness maatrix
        k_ele = k_ele + Bmat' * model.D * Bmat * detJ * weight(ll);

        % element nodal forces, body forces
        be = N * model.b(:, ee);   % interpolate body forces model.b to quadrature point, b(xi,eta) = N1(xi,eta) * b(:,ee) + ...
        f_ele = f_ele + N' * be * detJ * weight(ll);

    end % end of quadrature loop

        
end
