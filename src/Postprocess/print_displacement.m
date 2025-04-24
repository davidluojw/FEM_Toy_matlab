function print_displacement(model)
    % PRINT_DISPLACEMENT

    dis = zeros(model.neq, 1);
    for ii = 1:model.nnp
        for jj = 1:model.ndof
            pp = model.ndof*(ii - 1) + jj;
            dis(pp) = model.d(model.ID(jj, ii));
        end
    end
    
    % table title
    fprintf('\n                           Nodal displacement\n');
    fprintf('-------------------------------------------------------------------------------\n');
    fprintf('\tnode\tx\ty\t\tu_x\t\tu_y\t\tu_r\n');
    
    for i = 1:model.nnp
        x_coord = model.nodes(i, 1);
        y_coord = model.nodes(i, 2);
        u_x = dis(2*i - 1);  
        u_y = dis(2*i);   
        u_r = sqrt(u_x^2 + u_y^2);
        
        % print
        line = sprintf('\t%d\t%.2f\t%.2f\t%.15e\t%.15e\t%.15e', i, x_coord, y_coord, u_x, u_y, u_r);
        fprintf('%s\n', line);
    end
end