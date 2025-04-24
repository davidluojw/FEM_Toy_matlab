function model = solvedr(model)

    K = spalloc(model.neq, model.neq, 9*model.neq);

    K = model.K;

    f = model.f;

    % The "1" method
    for ii = 1:model.neq
        if model.flags(ii) == 2
            K(ii, :) = zeros(1, model.neq);
            K(:, ii) = zeros(model.neq, 1);
            K(ii, ii) = 1.0;
            f(ii) = model.e_bc(ii);
        else
            f(ii) = f(ii) - model.K(ii, :) * model.e_bc;
        end
    end

    model.d = K \ f;

    model.r = model.K * model.d - model.f;



end