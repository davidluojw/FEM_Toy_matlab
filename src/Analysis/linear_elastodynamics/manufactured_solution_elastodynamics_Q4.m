function model = manufactured_solution_elastodynamics_Q4(model)
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


    model.exact_ux    = @(x,y,t) sin(t)*x*(2-x)*y*(1-y); 
    model.exact_ux_x  = @(x,y,t) sin(t)*(2-2*x)*y*(1-y); 
    model.exact_ux_y  = @(x,y,t) sin(t)*x*(2-x)*(1-2*y); 
    model.exact_ux_xx = @(x,y,t) -sin(t)*2*y*(1-y)     ;  
    model.exact_ux_xy = @(x,y,t) sin(t)*(2-2*x)*(1-2*y);
    model.exact_ux_yx = @(x,y,t) sin(t)*(2-2*x)*(1-2*y);
    model.exact_ux_yy = @(x,y,t) -sin(t)*2*x*(2-x)     ;

    model.exact_uy    = @(x,y,t) cos(t)*x*(2-x)*y*(1-y);
    model.exact_uy_x  = @(x,y,t) cos(t)*(2-2*x)*y*(1-y);
    model.exact_uy_y  = @(x,y,t) cos(t)*x*(2-x)*(1-2*y);
    model.exact_uy_xx = @(x,y,t) -cos(t)*2*y*(1-y)     ;
    model.exact_uy_xy = @(x,y,t) cos(t)*(2-2*x)*(1-2*y);
    model.exact_uy_yx = @(x,y,t) cos(t)*(2-2*x)*(1-2*y);
    model.exact_uy_yy = @(x,y,t) -cos(t)*2*x*(2-x)     ;

    % velocity
    model.exact_vx = @(x,y,t) cos(t)*x*(2-x)*y*(1-y); 
    model.exact_vy = @(x,y,t) -sin(t)*x*(2-x)*y*(1-y);

    % acceleration
    model.exact_ax = @(x,y,t) -sin(t)*x*(2-x)*y*(1-y); 
    model.exact_ay = @(x,y,t) -cos(t)*x*(2-x)*y*(1-y);



    % strain
    model.exact_strainxx = @(x,y,t) model.exact_ux_x(x,y,t);
    model.exact_strainyy = @(x,y,t) model.exact_uy_y(x,y,t);
    model.exact_strainxy = @(x,y,t) 0.5 * (model.exact_ux_y(x,y,t) + model.exact_uy_x(x,y,t));

    % stress
    model.exact_stressxx = @(x,y,t) E / (1 - nu*nu) * ( model.exact_strainxx(x,y,t) + nu * model.exact_strainyy(x,y,t) );
    model.exact_stressyy = @(x,y,t) E / (1 - nu*nu) * ( nu * model.exact_strainxx(x,y,t) + model.exact_strainyy(x,y,t)  );
    model.exact_stressxy = @(x,y,t) E / (1+nu) * model.exact_strainxy(x,y,t);

    model.exact_stressxx_x = @(x,y,t) E / (1 - nu*nu) * (model.exact_ux_xx(x,y,t) + nu * model.exact_uy_yx(x,y,t) );
    model.exact_stressyy_y = @(x,y,t) E / (1 - nu*nu) * (nu * model.exact_ux_xy(x,y,t) + model.exact_uy_yy(x,y,t));
    model.exact_stressxy_x = @(x,y,t) E / (2*(1 + nu)) * (model.exact_uy_xx(x,y,t) + model.exact_ux_yx(x,y,t));
    model.exact_stressxy_y = @(x,y,t) E / (2*(1 + nu)) * (model.exact_uy_xy(x,y,t) + model.exact_ux_yy(x,y,t));

    % body force
    model.exact_bf_x = @(x,y,t) model.rho * model.exact_ax(x,y,t) - model.exact_stressxx_x(x,y,t) - model.exact_stressxy_y(x,y,t);
    model.exact_bf_y = @(x,y,t) model.rho * model.exact_ay(x,y,t) - model.exact_stressxy_x(x,y,t) - model.exact_stressyy_y(x,y,t);

    % traction force
    model.exact_t_x = @(x,y,t) model.exact_stressxy(x,y,t);
    model.exact_t_y = @(x,y,t) model.exact_stressyy(x,y,t);
    % model.exact_t_x = @(x,y) cos(atan(-y/x)) * model.exact_stressxx(x,y) + sin(atan(-y/x)) * model.exact_stressxy(x,y);
    % model.exact_t_y = @(x,y) cos(atan(-y/x)) * model.exact_stressxy(x,y) + sin(atan(-y/x)) * model.exact_stressyy(x,y);

    % Initial B.C 
    for ii = 1:model.nnp
        node_x = model.nodes(ii, 1);
        node_y = model.nodes(ii, 2);
        for jj = 1:model.ndof 
            dof_ids = model.ndof *(ii - 1) + jj;
            if jj == 1
                model.d_ini(dof_ids) = model.exact_ux(node_x, node_y, model.initial_t);
                model.v_ini(dof_ids) = model.exact_vx(node_x, node_y, model.initial_t);
            else
                model.d_ini(dof_ids) = model.exact_uy(node_x, node_y, model.initial_t);
                model.v_ini(dof_ids) = model.exact_vy(node_x, node_y, model.initial_t);
            end
        end
    end
    % store the disp
    model.disp(:, 1) = model.d_ini;


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

        % the deformed coordinates
        dis = zeros(model.neq, 1);
        for ii = 1:model.nnp
            node_x = model.nodes(ii, 1);
            node_y = model.nodes(ii, 2);

            dof_ids1 = ii * model.ndof - 1;
            dof_ids2 = ii * model.ndof;

            cur_t = model.time_intervals(tt);
            dis(dof_ids1) = model.exact_ux(node_x, node_y, cur_t) * model.fact;  % plot the result at time tt * dt
            dis(dof_ids2) = model.exact_uy(node_x, node_y, cur_t) * model.fact;
        end
        xnew = model.nodes(:,1) + dis(1:2:end);
        ynew = model.nodes(:,2) + dis(2:2:end);

        model.exact_disp(:, tt) = dis;


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
    movie(gcf, F, 1, 12); % play once, 15 frames per second

    % create video writing object
    v = VideoWriter('deformation.mp4', 'MPEG-4');
    open(v);

    for k = 1:numFrames
        writeVideo(v, F(k));
    end

    % close video file
    close(v);

end
