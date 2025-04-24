function [val_xi, val_eta] = Area_ShapeBasisN_Grad_2d(degree, a, xi, eta)

    psi = 1.0 - xi - eta;

    switch degree
        case 1
            if a == 1
                val_xi = 1.0;
                val_eta = 0.0;
            elseif a == 2
                val_xi = 0.0;
                val_eta = 1.0;
            elseif a == 3
                val_xi = -1.0;
                val_eta = -1.0;
            end
        
        case 2
            if a == 1
                val_xi = 4.0*xi - 1.0;
                val_eta = 0.0;
            elseif a == 2
                val_xi = 0.0;
                val_eta = 4.0*eta - 1.0;
            elseif a == 3 
                val_xi = -3.0 + 4.0*xi + 4.0*eta;
                val_eta = -3.0 + 4.0*xi + 4.0*eta;
            elseif a == 4
                val_xi = 4.0*eta;
                val_eta = 4.0*xi;
            elseif a == 5
                val_xi = -4.0*eta;
                val_eta = 4.0 - 4*xi -8.0*eta;
            elseif a == 6
                val_xi = 4.0 - 8.0*xi - 4.0*eta;
                val_eta = -4.0*xi;
            end
    
            
    otherwise
        error('Error: degree has to be 1, 2.');



end
