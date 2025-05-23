function plot_mesh_Q8(model)
    % PLOT_MESH plot the initial mesh   
    figure;
    hold on;
    
    %% plot the natural BC (red line)
    for i = 1:model.nbe
        node1 = model.n_bc(1, i);
        node2 = model.n_bc(2, i);
        node3 = model.n_bc(3, i);
        
        x1 = model.nodes(node1, 1);
        y1 = model.nodes(node1, 2);
        x2 = model.nodes(node2, 1);
        y2 = model.nodes(node2, 2);
        x3 = model.nodes(node3, 1);
        y3 = model.nodes(node3, 2);
        
        % plot the boundary
        plot([x1, x2, x3], [y1, y2, y3], 'Color','r', 'LineWidth',4);
    end
    
    %% plot the element mesh (blue line)
    for i = 1:model.nel
        node_ids = model.IEN(:,i);
        node_idsreo = [node_ids(1), node_ids(5), node_ids(2), node_ids(6), node_ids(3), node_ids(7), node_ids(4), node_ids(8)];  % reorder the index
        x_coords = model.nodes(node_idsreo, 1);
        y_coords = model.nodes(node_idsreo, 2);
        
        % plot the traingle
        XX = [x_coords; x_coords(1)];
        YY = [y_coords; y_coords(1)];
        plot(XX, YY, 'b');
        
        %% node number annotation
        for j = 1:model.nen
            text(x_coords(j), y_coords(j), num2str(node_idsreo(j)), 'FontSize',8, 'Color','k');
        end
    end
    
    title('Initial Structure');
    xlabel('X');
    ylabel('Y');
    axis equal;
    grid on;

    
    %% mesh parameters
    fprintf('\n=== Mesh Parameters ===\n');
    fprintf('No. of Elements:  %d\n', model.nel);
    fprintf('No. of Nodes:     %d\n', model.nnp);
    fprintf('No. of Equations: %d\n\n', model.neq);

    hold off;

end