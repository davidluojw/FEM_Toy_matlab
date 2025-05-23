function plot_displacement_T3(model)
    % PLOT THE DEFORMED DISPLACEMENT 
    dis = zeros(model.neq, 1);
    for ii = 1:model.nnp
        for jj = 1:model.ndof
            pp = model.ndof*(ii - 1) + jj;
            dis(pp) = model.d(model.ID(jj, ii)) * model.fact;
        end
    end

    
    %% the deformed coordinates
    xnew = model.nodes(:, 1) + dis(1:2:end);  
    ynew = model.nodes(:, 2) + dis(2:2:end);  
    
    figure;
    hold on;
    for i = 1:model.nel
        % initial structure (blue line)
        node_ids = model.IEN(:,i);
        X_orig = model.nodes(node_ids([1:end,1]), 1);  
        Y_orig = model.nodes(node_ids([1:end,1]), 2);

        if i == 1
            h_initial = plot(X_orig, Y_orig, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Initial');
        else
            plot(X_orig, Y_orig, 'b-', 'LineWidth', 1.5);
        end
    end

    %% the deformed structure (black line)
    for i = 1:model.nel
        node_ids = model.IEN(:,i);
        X_def = xnew(node_ids([1:end,1]));
        Y_def = ynew(node_ids([1:end,1]));
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

     %% the displacement  contour
     plot_displacement_contour_T3(model, 1);  % u_x
     plot_displacement_contour_T3(model, 2);  % u_y
     plot_displacement_contour_T3(model, 3);  % u_abs
 

end