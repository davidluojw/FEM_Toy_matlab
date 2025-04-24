function [xi, eta, w] = gauss_2d(N1, N2)
    % N1: the number of quadrature points for xi-coor
    % N2: the number of quadrature points for eta-coor

    % preallocation
    xi = zeros(N1*N2,1);
    eta = xi;
    w = xi;
    
    % generate 1D rule
    [x1, w1] = gauss(N1, -1, 1);
    
    [x2, w2] = gauss(N2, -1, 1);
    
    for ii = 1 : N1
        for jj = 1 : N2
        xi( (jj-1)*N1 + ii )  = x1(ii);
        eta( (jj-1)*N1 + ii ) = x2(jj);
        w( (jj-1)*N1 + ii)    = w1(ii) * w2(jj);
        end
    end
    
    % EOF