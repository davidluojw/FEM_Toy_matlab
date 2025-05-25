function [val_xi, val_eta] = T6_ShapeBasisN_Grad_2d(a, xi, eta)
    % a: the number of element node in xi-eta coordinate Q8
    % xi: the quadrature point value of xi-coor
    % eta: the quadrature point value of eta-coor

    switch a
        case 1 
            val_xi = 4*xi - 1;
            val_eta = 0;
        
        case 2
            val_xi = 0;
            val_eta = 4*eta - 1;

        case 3
            val_xi = 4*xi + 4*eta - 3;
            val_eta = 4*xi + 4*eta - 3;

        case 4
            val_xi = 4*eta;
            val_eta = 4*xi;

        case 5
            val_xi = -4*eta;
            val_eta = 4*(1-xi-2*eta);
        
        case 6
            val_xi = 4*(1-eta-2*xi);
            val_eta = -4*xi;


        otherwise
            error('Error: a has to be 1, 2, ...  or 6.');

end