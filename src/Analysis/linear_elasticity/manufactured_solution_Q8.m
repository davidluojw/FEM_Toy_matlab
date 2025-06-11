function model = manufactured_solution_Q8(model)
    E = model.E;
    nu = model.nu;

    G    = @(x, y) sin((x + y) * 2 * pi);
    G_x  = @(x, y) 2 * pi * cos((x + y) * 2 * pi);
    G_y  = @(x, y) 2 * pi * cos((x + y) * 2 * pi);
    G_xx = @(x, y) -4 * pi * pi * sin((x + y) * 2 * pi);
    G_yy = @(x, y) -4 * pi * pi * sin((x + y) * 2 * pi);
    G_xy = @(x, y) -4 * pi * pi * sin((x + y) * 2 * pi);
    G_yx = @(x, y) -4 * pi * pi * sin((x + y) * 2 * pi);

    % model.exact_ux    = @(x,y) x*(2-x)*y*(1-y) + 0.1*G(x, y);
    % model.exact_ux_x  = @(x,y) (2-2*x)*y*(1-y) + 0.1*G_x(x, y); 
    % model.exact_ux_y  = @(x,y) x*(2-x)*(1-2*y) + 0.1*G_y(x, y);
    % model.exact_ux_xx = @(x,y) -2*y*(1-y)      + 0.1*G_xx(x, y); 
    % model.exact_ux_xy = @(x,y) (2-2*x)*(1-2*y) + 0.1*G_xy(x, y); 
    % model.exact_ux_yx = @(x,y) (2-2*x)*(1-2*y) + 0.1*G_yx(x, y);
    % model.exact_ux_yy = @(x,y) -2*x*(2-x)      + 0.1*G_yy(x, y);
    % 
    % model.exact_uy    = @(x,y) x*(2-x)*y*(1-y) + 0.1*G(x, y);
    % model.exact_uy_x  = @(x,y) (2-2*x)*y*(1-y) + 0.1*G_x(x, y);
    % model.exact_uy_y  = @(x,y) x*(2-x)*(1-2*y) + 0.1*G_y(x, y);
    % model.exact_uy_xx = @(x,y) -2*y*(1-y)      + 0.1*G_xx(x, y);
    % model.exact_uy_xy = @(x,y) (2-2*x)*(1-2*y) + 0.1*G_xy(x, y);
    % model.exact_uy_yx = @(x,y) (2-2*x)*(1-2*y) + 0.1*G_yx(x, y);
    % model.exact_uy_yy = @(x,y) -2*x*(2-x)      + 0.1*G_yy(x, y);


    model.exact_ux    = @(x,y) x*(2-x)*y*(1-y); 
    model.exact_ux_x  = @(x,y) (2-2*x)*y*(1-y); 
    model.exact_ux_y  = @(x,y) x*(2-x)*(1-2*y); 
    model.exact_ux_xx = @(x,y) -2*y*(1-y)     ;  
    model.exact_ux_xy = @(x,y) (2-2*x)*(1-2*y);
    model.exact_ux_yx = @(x,y) (2-2*x)*(1-2*y);
    model.exact_ux_yy = @(x,y) -2*x*(2-x)     ;

    model.exact_uy    = @(x,y) x*(2-x)*y*(1-y);
    model.exact_uy_x  = @(x,y) (2-2*x)*y*(1-y);
    model.exact_uy_y  = @(x,y) x*(2-x)*(1-2*y);
    model.exact_uy_xx = @(x,y) -2*y*(1-y)     ;
    model.exact_uy_xy = @(x,y) (2-2*x)*(1-2*y);
    model.exact_uy_yx = @(x,y) (2-2*x)*(1-2*y);
    model.exact_uy_yy = @(x,y) -2*x*(2-x)     ;



    % strain
    model.exact_strainxx = @(x,y) model.exact_ux_x(x,y);
    model.exact_strainyy = @(x,y) model.exact_uy_y(x,y);
    model.exact_strainxy = @(x,y) 0.5 * (model.exact_ux_y(x,y) + model.exact_uy_x(x,y));

    % stress
    model.exact_stressxx = @(x,y) E / (1 - nu*nu) * ( model.exact_strainxx(x,y) + nu * model.exact_strainyy(x,y) );
    model.exact_stressyy = @(x,y) E / (1 - nu*nu) * ( nu * model.exact_strainxx(x,y) + model.exact_strainyy(x,y)  );
    model.exact_stressxy = @(x,y) E / (1+nu)* model.exact_strainxy(x,y);

    model.exact_stressxx_x = @(x,y) E / (1 - nu*nu) * (model.exact_ux_xx(x,y) + nu * model.exact_uy_yx(x,y) );
    model.exact_stressyy_y = @(x,y) E / (1 - nu*nu) * (nu * model.exact_ux_xy(x,y) + model.exact_uy_yy(x,y));
    model.exact_stressxy_x = @(x,y) E / (2*(1 + nu)) * (model.exact_uy_xx(x,y) + model.exact_ux_yx(x,y));
    model.exact_stressxy_y = @(x,y) E / (2*(1 + nu)) * (model.exact_uy_xy(x,y) + model.exact_ux_yy(x,y));

    % body force
    model.exact_bf_x = @(x,y) -model.exact_stressxx_x(x,y) - model.exact_stressxy_y(x,y);
    model.exact_bf_y = @(x,y) -model.exact_stressxy_x(x,y) - model.exact_stressyy_y(x,y);

    % traction force
    model.exact_t_x = @(x,y) model.exact_stressxy(x,y);
    model.exact_t_y = @(x,y) model.exact_stressyy(x,y);
    % model.exact_t_x = @(x,y) cos(atan(-y/x)) * model.exact_stressxx(x,y) + sin(atan(-y/x)) * model.exact_stressxy(x,y);
    % model.exact_t_y = @(x,y) cos(atan(-y/x)) * model.exact_stressxy(x,y) + sin(atan(-y/x)) * model.exact_stressyy(x,y);

    % Dirichlet B.C. 
    for ii = 1:model.nnp
        node_x = model.nodes(ii, 1);
        node_y = model.nodes(ii, 2);
        for jj = 1:model.ndof 
            dof_ids = model.ndof *(ii - 1) + jj;
            if model.flags(dof_ids) == 2 &&  jj == 1
                model.e_bc(dof_ids) = model.exact_ux(node_x, node_y);
            elseif model.flags(dof_ids) == 2 &&  jj == 2
                model.e_bc(dof_ids) = model.exact_uy(node_x, node_y);
            end
        end
    end

    % Neumann B.C
    % model.nbc_nodes = [21, 22, 23, 24; 22, 23, 24, 25];
    % model.nbc_nodes = [1, 2, 3, 4; 2, 3, 4, 5];


    % PLOT THE DEFORMED DISPLACEMENT (exact solution)
    dis = zeros(model.neq, 1);
    for ii = 1:model.nnp
        node_x = model.nodes(ii, 1);
        node_y = model.nodes(ii, 2);

        dof_ids1 = model.ID(1, ii);
        dof_ids2 = model.ID(2, ii);
        
        dis(dof_ids1) = model.exact_ux(node_x, node_y) * model.fact;
        dis(dof_ids2) = model.exact_uy(node_x, node_y) * model.fact;

    end

    
    %% the deformed coordinates
    xnew = model.nodes(:, 1) + dis(1:2:end);  
    ynew = model.nodes(:, 2) + dis(2:2:end);  
    
    figure;
    hold on;
    for i = 1:model.nel
        % initial structure (blue line)
        node_ids = model.IEN(:,i);
        node_idsreo = [node_ids(1), node_ids(5), node_ids(2), node_ids(6), node_ids(3), node_ids(7), node_ids(4), node_ids(8)];  % reorder the index
        X_orig = model.nodes(node_idsreo([1:end,1]), 1);  
        Y_orig = model.nodes(node_idsreo([1:end,1]), 2);

        if i == 1
            h_initial = plot(X_orig, Y_orig, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Initial');
        else
            plot(X_orig, Y_orig, 'b-', 'LineWidth', 1.5);
        end
    end

    %% the deformed structure (black line)
    for i = 1:model.nel
        node_ids = model.IEN(:,i);
        node_idsreo = [node_ids(1), node_ids(5), node_ids(2), node_ids(6), node_ids(3), node_ids(7), node_ids(4), node_ids(8)];  % reorder the index
        X_def = xnew(node_idsreo([1:end,1]));
        Y_def = ynew(node_idsreo([1:end,1]));
        if i == 1
            h_deformed = plot(X_def, Y_def, 'k--', 'LineWidth', 2, 'DisplayName', 'Deformed');
        else
            plot(X_def, Y_def, 'k--', 'LineWidth', 2);
        end
    end
    
    title('Initial and Deformed Structure');
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    legend([h_initial, h_deformed], 'Location', 'best'); 
    axis equal;
    grid on;
    hold off;

    

end
