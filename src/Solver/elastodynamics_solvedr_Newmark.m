function model = elastodynamics_solvedr_Newmark(model, tt)

    LEFT = spalloc(model.neq, model.neq, (model.ndof * model.nen + 1)*model.neq);

    K = model.K;

    M = model.M;

    f = model.f;

    % Dirichlet B.C. 
    cur_t = model.time_intervals(tt);
    for ii = 1:model.nnp
        node_x = model.nodes(ii, 1);
        node_y = model.nodes(ii, 2);
        for jj = 1:model.ndof 
            dof_ids = model.ndof *(ii - 1) + jj;
            if model.flags(dof_ids) == 2 &&  jj == 1
                model.e_bc(dof_ids) = model.exact_ux(node_x, node_y, cur_t);
            elseif model.flags(dof_ids) == 2 &&  jj == 2
                model.e_bc(dof_ids) = model.exact_uy(node_x, node_y, cur_t);
            end
        end
    end

    if tt <= 2

        % acceleration 
        model.d_ddot = model.M \ (model.f - model.K * model.d_ini);


        % effective stiffness matrix
        K = 1 / (model.beta * model.dt^2) * model.M + model.K;
        LEFT = K;

        % effective force vector
        RIGHT = model.f + 1/(model.beta * model.dt^2) * model.M * model.d_ini + 1/(model.beta * model.dt) * model.M * model.v_ini + (1/(2*model.beta) - 1) * model.M * model.d_ddot; 

        model.d_dot = model.v_ini;

    else
        % effective stiffness matrix
        K = 1 / (model.beta * model.dt^2) * model.M + model.K;
        LEFT = K;

        % effective force vector
        RIGHT = model.f +  1/(model.beta * model.dt^2) * model.M * model.disp(:, tt-1) + 1/(model.beta * model.dt) * model.M * model.d_dot + (1/(2*model.beta) - 1) * model.M * model.d_ddot; 

    end

    % The "1" method
    for ii = 1:model.neq
        if model.flags(ii) == 2
            LEFT(ii, :) = zeros(1, model.neq);
            LEFT(:, ii) = zeros(model.neq, 1);
            LEFT(ii, ii) = 1.0;
            RIGHT(ii) = model.e_bc(ii);
        else
            RIGHT(ii) = RIGHT(ii) - K(ii, :) * model.e_bc;
        end
    end

    % disp at t = 1
    model.d = LEFT \ RIGHT;

    model.r = model.K * model.d - model.f;

    % store the solution 
    model.disp(:, tt) = model.d; 

    % update the acceleration and velocity
    model.d_dot = model.gamma/(model.beta * model.dt) * (model.disp(:, tt) - model.disp(:, tt-1))  + (1 - model.gamma/model.beta) * model.d_dot + (1 - model.gamma / (2*model.beta)) * model.dt * model.d_ddot;
    model.d_ddot = 1/(model.beta * model.dt^2) * (model.disp(:, tt) - model.disp(:, tt-1)) - 1/(model.beta * model.dt) * model.d_dot - (1/(2*model.beta) - 1) * model.d_ddot;



end