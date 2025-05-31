function play_displacement_Q4(model) 

    % PLOT THE DEFORMED DISPLACEMENT (exact solution)
    % Initialize animation
    figure;
    hold on;
    axis equal;
    grid on;
    xlabel('X'); ylabel('Y'); 
    title('Structural Deformation Animation');

    numFrames = model.nts;   % number of frames
    F(numFrames) = struct('cdata', [], 'colormap', []);

    model.exact_disp = zeros(model.neq, model.nts);

    %% 主循环
    for tt = 1:numFrames
        cla(gca); % clear current axis [1,2](@ref)
        hold on;  % re-enable graphic preservation

        dis = zeros(model.neq, 1);
        for ii = 1:model.nnp
            for jj = 1:model.ndof
                pp = model.ndof*(ii - 1) + jj;
                dis(pp) = model.disp(model.ID(jj, ii), tt) * model.fact;
            end
        end

        xnew = model.nodes(:, 1) + dis(1:2:end);  
        ynew = model.nodes(:, 2) + dis(2:2:end);  


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

        %  the deformed structure (black line)
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

        % update the graph
        legend([h_initial, h_deformed], 'Location', 'best'); 
        drawnow limitrate;
        pause(0.005);
        F(tt) = getframe(gcf);

    end

    % play animation
    movie(gcf, F, 1, 120); % play once, 15 frames per second

    % create video writing object
    v = VideoWriter('play_deformation_Q4.mp4', 'MPEG-4');
    v.FrameRate = 120;
    open(v);

    for k = 1:numFrames
        writeVideo(v, F(k));
    end

    % close video file
    close(v);


end