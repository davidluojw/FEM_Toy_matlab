function val = Q8_ShapeBasisN_2d(a, xi, eta)
    % a: the number of element node in xi-eta coordinate Q8
    % xi: the quadrature point value of xi-coor
    % eta: the quadrature point value of eta-coor

    switch a
        case 1 
            val = 0.25*(1-xi)*(eta-1)*(xi+eta+1);
        
        case 2
            val = 0.25*(1+xi)*(eta-1)*(eta-xi+1);

        case 3
            val = 0.25*(1+xi)*(1+eta)*(xi+eta-1);

        case 4
            val = 0.25*(xi-1)*(eta+1)*(xi-eta+1);

        case 5
            val = 0.5*(1-eta)*(1-xi*xi);
        
        case 6
            val = 0.5*(1+xi)*(1-eta*eta);

        case 7
            val = 0.5*(1+eta)*(1-xi*xi);

        case 8
            val = 0.5*(1-xi)*(1-eta*eta);


        otherwise
            error('Error: a has to be 1, 2, ...  or 8.');

end