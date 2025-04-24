function plot_displacement_contour_T3(model, type)

    %% the displacement u_x contour

    dis = zeros(model.neq, 1);
    for ii = 1:model.nnp
        for jj = 1:model.ndof
            pp = model.ndof*(ii - 1) + jj;
            dis(pp) = model.d(model.ID(jj, ii)) * model.fact;
        end
    end


    figure;
    hold on;
    
    % color range
    % vmin = -150;
    % vmax = 200;
    plot_option();
    
    % loop over elements
    for i = 1:model.nel
        node_ids = model.IEN(:,i);
        x_nodes = model.nodes(node_ids, 1);
        y_nodes = model.nodes(node_ids, 2);
        
        % Quadrilateral nodal coordinates
        XX = [x_nodes(1), x_nodes(2), x_nodes(3), x_nodes(1)];
        YY = [y_nodes(1), y_nodes(2), y_nodes(3), y_nodes(1)];
        
        % nodal displacement
        u_x = dis(model.ID(1, node_ids));
        u_y = dis(model.ID(2, node_ids));
        uu = sqrt(u_x .* u_x + u_y .* u_y);
        
        % stress matrix
        if type == 1
            dd = [u_x(1), u_x(2), u_x(3), u_x(1)];
        elseif type == 2
            dd = [u_y(1), u_y(2), u_y(3), u_y(1)];
        elseif type == 3
            dd = [uu(1), uu(2), uu(3), uu(1)];
        end

        patch(XX,YY,dd);
        hold on;
    end

    % colormap(jet);
    h_cb  = colorbar('Location', 'eastoutside');
    h_cb.Label.String = 'Displacement ';
    
    % plot mesh
    for i = 1:model.nel
        node_ids = model.IEN(:,i);
        X_contour = model.nodes(node_ids, 1); 
        Y_contour = model.nodes(node_ids, 2);
        plot(X_contour, Y_contour, 'k-', 'LineWidth', 0.5);
    end
    
    if type == 1
        title('Displacement U_x Contours');
    elseif type == 2
        title('Displacement U_y Contours');
    elseif type == 3
        title('Displacement U_{abs} Contours');
    end
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    axis equal;
    
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 6]);
    print('elasticity-sxx', '-dpdf', '-r300');
    hold off;


end