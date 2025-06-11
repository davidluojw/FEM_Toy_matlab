clear; clc; close all

currentDir = pwd;

addpath(genpath(currentDir));

file_in = fullfile(currentDir, 'input', 'input_2d_Q4_ex4_5.dat'); 

fem_data = read_fem_Q4_dat(file_in);

model = create_2d_Q4_model(fem_data);

model = setup_ID_LM(model);

plot_mesh(model);

% Use the manufactured solution
model = elastodynamics_manufactured_solution_Q4(model);   % Update model.e_bc and model.n_bc


for tt = 1:model.nts-1

    model.f = zeros(model.neq, 1);

    for ee = 1:model.nel
        model = elastodynamics_elemQuad4_2d(model, ee, tt);
        model = elastodynamics_assembly(model, ee, tt);
    end
    
    model = elastodynamics_concentrate_traction_Q4(model, tt+1); 
    
    % model = elastodynamics_solvedr_Central(model, tt+1); 
    model = elastodynamics_solvedr_Newmark(model, tt+1); 
    % model = elastodynamics_solvedr_GenAlpha(model, tt+1); 

    print_displacement(model);

    % plot_displacement_Q4(model);
    % model = get_qdpt_stress_Q4(model);
    % model = get_nodal_stress_Q4(model);
    % cart2polar_stress(model);
    % plot_stress_Q4(model);
    
end

play_displacement_Q4(model);
model = elastodynamics_error_analysis_Q4(model);



