function val = T6_ShapeBasisN_2d(a, xi, eta, zeta)
    % a: the number of element node in xi-eta coordinate Q8
    % xi: the quadrature point value of xi-coor
    % eta: the quadrature point value of eta-coor

    switch a
        case 1 
            val = xi*(2*xi - 1);
        
        case 2
            val = eta*(2*eta - 1);

        case 3
            val = zeta*(2*zeta - 1);

        case 4
            val = 4*xi*eta;

        case 5
            val = 4*eta*zeta;
        
        case 6
            val = 4*xi*zeta;


        otherwise
            error('Error: a has to be 1, 2, ...  or 6.');

end