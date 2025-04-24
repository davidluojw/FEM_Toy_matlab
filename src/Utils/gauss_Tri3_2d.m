function [xi, eta, zeta, weight] = gauss_Tri3_2d(N)
    % the number of quadrature ppints

    switch N
        case 1
            xi = [1/3];
            eta = [1/3];
            zeta = [1/3];
            weight = [1/2];
        case 2
            xi = [2/3, 1/6, 1/6];
            eta = [1/6, 2/3, 1/6];
            zeta = [1/6, 1/6, 2/3];
            weight = [1/6, 1/6, 1/6];

        case 3
            xi = [1/2, 1/2, 0];
            eta = [1/2, 0, 1/2];
            zeta = [0, 1/2, 1/2];
            weight = [1/6, 1/6, 1/6];

        case 4
            xi = [1/3, 0.6, 0.2, 0.2];
            eta = [1/3, 0.2, 0.6, 0.2];
            zeta = [1/3, 0.2, 0.2, 0.6];
            weight = [-27/96, 25/96, 25/96, 25/96];
        
        
    otherwise
        error('Error: N has to be 1, 2, 3, or 4.');



end
