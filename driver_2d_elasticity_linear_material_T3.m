clear; clc

currentDir = pwd;

addpath(genpath(currentDir));

file_in = fullfile(currentDir, 'input', 'input_2d_T3_.dat'); 

fem_data = read_fem_T3_dat(file_in);

model = create_2d_T3_model(fem_data);

model = setup_ID_LM(model);

for ee = 1:model.nel
    [k_ele, f_ele] = elasticity_elemTri3_2d(model, ee);
    model = assembly(model, k_ele, f_ele, ee);
end

model = concentrate_traction(model);

model = solvedr(model);

plot_mesh(model);

print_displacement(model);

plot_displacement_T3(model);

model = get_qdpt_stress_T3(model);
model = get_nodal_stress_T3(model);


cart2polar_stress(model);

model = plot_stress_T3(model);