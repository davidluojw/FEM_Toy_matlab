function plot_stress_contour_Q4(model, type)

    %% von mmises stress
    figure;
    hold on;
    
    % color range
    % vmin = 0;
    % vmax = 250;
    plot_option();
    
    for i = 1:model.nel
        node_ids = model.IEN(:,i);
        x_nodes = model.nodes(node_ids, 1);
        y_nodes = model.nodes(node_ids, 2);
        
        XX = [x_nodes(1), x_nodes(2), x_nodes(3), x_nodes(4), x_nodes(1)];
        YY = [y_nodes(1), y_nodes(2), y_nodes(3), y_nodes(4), y_nodes(1)];
        
        % each componenet of stress
        sxx = model.stress_nodal(1, node_ids);
        syy = model.stress_nodal(2, node_ids);
        sxy = model.stress_nodal(3, node_ids);
        
        % von mises stress
        S1 = 0.5*(sxx + syy) + sqrt((0.5*(sxx - syy)).^2 + sxy.^2);
        S2 = 0.5*(sxx + syy) - sqrt((0.5*(sxx - syy)).^2 + sxy.^2);
        mises = sqrt(S1.^2 + S2.^2 - S1.*S2);
        
        if type == 1
            dd = [sxx(1), sxx(2), sxx(3), sxx(4), sxx(1)];
        elseif type == 2
            dd = [syy(1), syy(2), syy(3), syy(4), syy(1)];
        elseif type == 3
            dd = [sxy(1), sxy(2), sxy(3), sxy(4), sxy(1)];
        elseif type == 4
            dd = [mises(1), mises(2), mises(3), mises(4), mises(1)];
        end
        patch(XX,YY,dd);
        hold on;
    end
    
    % colormap(jet);
    % caxis([vmin, vmax]);
    h_cb  = colorbar('Location', 'eastoutside');
    h_cb.Label.String = 'Stress ';
    
    % plot mesh 
    for i = 1:model.nel
        node_ids = model.IEN(:,i);
        X_contour = model.nodes(node_ids, 1);
        Y_contour = model.nodes(node_ids, 2);
        plot(X_contour, Y_contour, 'k-', 'LineWidth', 0.5);
    end
    if type == 1
        title('\sigma_{xx} Stress Contourss');
    elseif type == 2
        title('\sigma_{yy} Stress Contourss');
    elseif type == 3
        title('\sigma_{xy} Stress Contourss');
    elseif type == 4
        title('Von Mises Stress Contours');
    end
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    axis equal;

    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 6]);
    print('elasticity-mises', '-dpdf', '-r300');
    hold off;


end
