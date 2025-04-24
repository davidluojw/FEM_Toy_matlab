function val = Area_ShapeBasisN_2d(degree, a, xi, eta, zeta)

    switch degree
        case 1
            if a == 1
                val = xi;
            elseif a == 2
                val = eta;
            elseif a == 3
                val = zeta;
            end
        
        case 2
            if a == 1
                val = xi*(2.0*xi - 1.0);
            elseif a == 2
                val = eta*(2.0*eta - 1.0);
            elseif a == 3 
                val = zeta*(2.0*zeta - 1.0);
            elseif a == 4
                val = 4.0*xi*eta;
            elseif a == 5
                val = 4.0*eta*zeta;
            elseif a == 6
                val = 4.0*xi*zeta;
            end
    

    otherwise
        error('Error: degree has to be 1, 2.');
    



end
