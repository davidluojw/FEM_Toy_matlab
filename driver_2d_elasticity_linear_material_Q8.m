clear; clc

currentDir = pwd;

addpath(genpath(currentDir));

file_in = fullfile(currentDir, 'input', 'input_2d_Q8_ex4_5.dat'); 

fem_data = read_fem_Q8_dat(file_in);

model = create_2d_Q8_model(fem_data);

model = setup_ID_LM(model);

% Use the manufactured solution
model = manufactured_solution_Q8(model);   % Update model.e_bc and model.n_bc

for ee = 1:model.nel
    [k_ele, f_ele] = elasticity_elemQuad8_2d(model, ee);
    model = assembly(model, k_ele, f_ele, ee);
end

model = concentrate_traction_Q8(model); 



model = solvedr(model); 

print_displacement(model);

plot_mesh_Q8(model);

plot_displacement_Q8(model);

model = error_analysis_Q8(model);

model = get_qdpt_stress_Q8(model);
model = get_nodal_stress_Q8(model);

cart2polar_stress(model);

plot_stress_Q8(model);

