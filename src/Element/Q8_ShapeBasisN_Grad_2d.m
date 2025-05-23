function [val_xi, val_eta] = Q8_ShapeBasisN_Grad_2d(a, xi, eta)
    % a: the number of element node in xi-eta coordinate Q8
    % xi: the quadrature point value of xi-coor
    % eta: the quadrature point value of eta-coor

    switch a
        case 1 
            val_xi = 0.25*(1-eta)*(2*xi+eta);
            val_eta = 0.25*(1-xi)*(xi+2*eta);
        
        case 2
            val_xi = 0.25*(1-eta)*(2*xi-eta);
            val_eta = 0.25*(1+xi)*(2*eta-xi);

        case 3
            val_xi = 0.25*(1+eta)*(2*xi+eta);
            val_eta = 0.25*(1+xi)*(2*eta+xi);

        case 4
            val_xi = 0.25*(1+eta)*(2*xi-eta);
            val_eta = 0.25*(1-xi)*(2*eta-xi);

        case 5
            val_xi = xi*(eta-1);
            val_eta = 0.5*(xi*xi-1);
        
        case 6
            val_xi = 0.5*(1-eta*eta);
            val_eta = -eta*(xi+1);

        case 7
            val_xi = -xi*(eta+1);
            val_eta = 0.5*(1-xi*xi);

        case 8
            val_xi = 0.5*(eta*eta-1);
            val_eta = (xi-1)*eta;


        otherwise
            error('Error: a has to be 1, 2, ...  or 8.');

end