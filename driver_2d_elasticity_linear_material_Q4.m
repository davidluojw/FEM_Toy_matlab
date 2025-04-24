clear; clc

currentDir = pwd;

addpath(genpath(currentDir));

file_in = fullfile(currentDir, 'input', 'input_2d_Q4_Cyn_6_12.dat'); 

fem_data = read_fem_Q4_dat(file_in);

model = create_2d_Q4_model(fem_data);

model = setup_ID_LM(model);

for ee = 1:model.nel
    [k_ele, f_ele] = elasticity_elemQuad4_2d(model, ee);
    model = assembly(model, k_ele, f_ele, ee);
end

model = concentrate_traction(model);
model = solvedr(model); 



print_displacement(model);

model = get_qdpt_stress_Q4(model);
model = get_nodal_stress_Q4(model);

cart2polar_stress(model);

plot_mesh(model);

plot_displacement_Q4(model);

model = plot_stress_Q4(model);
