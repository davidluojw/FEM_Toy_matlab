clear; clc

currentDir = pwd;

addpath(genpath(currentDir));

file_in = fullfile(currentDir, 'input', 'input_2d_T6_ex4_1.dat'); 

fem_data = read_fem_T6_dat(file_in);

model = create_2d_T6_model(fem_data);

model = setup_ID_LM(model);

% Use the manufactured solution
model = manufactured_solution_T6(model);   % Update model.e_bc and model.n_bc

for ee = 1:model.nel
    [k_ele, f_ele] = elasticity_elemTri6_2d(model, ee);
    model = assembly(model, k_ele, f_ele, ee);
end

model = concentrate_traction_T6(model); 

model = solvedr(model); 

print_displacement(model);

plot_mesh_T6(model);

plot_displacement_T6(model);

model = error_analysis_T6(model);

model = get_qdpt_stress_T6(model);
model = get_nodal_stress_T6(model);

cart2polar_stress(model);

plot_stress_T6(model);


