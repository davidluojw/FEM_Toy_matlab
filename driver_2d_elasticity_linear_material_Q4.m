clear; clc

currentDir = pwd;

addpath(genpath(currentDir));

file_in = fullfile(currentDir, 'input', 'input_2d_Q4_ex4_5.dat'); 
% file_in = fullfile(currentDir, 'input/genInput', 'plate_mesh_8x8.dat'); 

fem_data = read_fem_Q4_dat(file_in);

model = create_2d_Q4_model(fem_data);

model = setup_ID_LM(model);

plot_mesh(model);

% Use the manufactured solution
model = manufactured_solution_Q4(model);   % Update model.e_bc and model.n_bc

for ee = 1:model.nel
    [k_ele, f_ele] = elasticity_elemQuad4_2d(model, ee);
    model = assembly(model, k_ele, f_ele, ee);
end

model = concentrate_traction_Q4(model); 

model = solvedr(model); 

model = error_analysis_Q4(model);

print_displacement(model);

plot_displacement_Q4(model);

model = get_qdpt_stress_Q4(model);
model = get_nodal_stress_Q4(model);

cart2polar_stress(model);

plot_stress_Q4(model);
