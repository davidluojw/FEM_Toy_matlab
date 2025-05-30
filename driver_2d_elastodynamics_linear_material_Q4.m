clear; clc

currentDir = pwd;

addpath(genpath(currentDir));

file_in = fullfile(currentDir, 'input', 'input_2d_Q4_ex4_1.dat'); 

fem_data = read_fem_Q4_dat(file_in);

model = create_2d_Q4_model(fem_data);

model = setup_ID_LM(model);

plot_mesh(model);

% Use the manufactured solution
model = manufactured_solution_elastodynamics_Q4(model);   % Update model.e_bc and model.n_bc



for tt = 1:model.nts-1

    model.f = zeros(model.neq, 1);

    for ee = 1:model.nel
        model = elastodynamics_elemQuad4_2d(model, ee, tt);
        model = elastodynamics_assembly(model, ee, tt);
    end
    
    model = elastodynamics_concentrate_traction_Q4(model, tt+1); 
    
    model = elastodynamics_solvedr(model, tt+1); 
    
end



model = error_analysis_Q4(model);
    
print_displacement(model);

plot_displacement_Q4(model);

model = get_qdpt_stress_Q4(model);
model = get_nodal_stress_Q4(model);

cart2polar_stress(model);

plot_stress_Q4(model);

