FEM_Toy_matlab
--------------

# Program Structure
```plaintext
FEM/
├── Mesh_Cylinder.f90
├── README.md
├── driver_2d_elasticity_linear_material_Q4.m
├── driver_2d_elasticity_linear_material_Q8.m
├── driver_2d_elasticity_linear_material_T3.m
├── driver_2d_elasticity_linear_material_T6.m
├── elasticity-mises.pdf
├── elasticity-sxx.pdf
├── input
└── src
    ├── Analysis
    │   ├── convergence_analysis_Q4.m
    │   ├── convergence_analysis_Q8.m
    │   ├── error_analysis_Q4.m
    │   ├── error_analysis_Q8.m
    │   ├── error_analysis_T3.m
    │   ├── error_analysis_T6.m
    │   ├── manufactured_solution_Q4.m
    │   ├── manufactured_solution_Q8.m
    │   ├── manufactured_solution_T3.m
    │   └── manufactured_solution_T6.m
    ├── BoundaryCondition
    │   ├── concentrate_traction_Q4.m
    │   ├── concentrate_traction_Q8.m
    │   ├── concentrate_traction_T3.m
    │   └── concentrate_traction_T6.m
    ├── Element
    │   ├── Area_ShapeBasisN_2d.m
    │   ├── Area_ShapeBasisN_Grad_2d.m
    │   ├── Poly_ShapeBasisN_1d.m
    │   ├── Q8_ShapeBasisN_2d.m
    │   ├── Q8_ShapeBasisN_Grad_2d.m
    │   ├── Quad_ShapeBasisN_2d.m
    │   ├── Quad_ShapeBasisN_Grad_2d.m
    │   ├── T6_ShapeBasisN_2d.m
    │   └── T6_ShapeBasisN_Grad_2d.m
    ├── Model
    │   ├── elasticity_elemQuad4_2d.m
    │   ├── elasticity_elemQuad8_2d.m
    │   ├── elasticity_elemTri3_2d.m
    │   └── elasticity_elemTri6_2d.m
    ├── Postprocess
    │   ├── cart2polar_stress.m
    │   ├── get_nodal_stress_Q4.m
    │   ├── get_nodal_stress_Q8.m
    │   ├── get_nodal_stress_T3.m
    │   ├── get_nodal_stress_T6.m
    │   ├── get_qdpt_stress_Q4.m
    │   ├── get_qdpt_stress_Q8.m
    │   ├── get_qdpt_stress_T3.m
    │   ├── get_qdpt_stress_T6.m
    │   ├── plot_displacement_Q4.m
    │   ├── plot_displacement_Q8.m
    │   ├── plot_displacement_T3.m
    │   ├── plot_displacement_T6.m
    │   ├── plot_displacement_contour_Q4.m
    │   ├── plot_displacement_contour_Q8.m
    │   ├── plot_displacement_contour_T3.m
    │   ├── plot_displacement_contour_T6.m
    │   ├── plot_mesh.m
    │   ├── plot_mesh_Q8.m
    │   ├── plot_mesh_T6.m
    │   ├── plot_stress_Q4.m
    │   ├── plot_stress_Q8.m
    │   ├── plot_stress_T3.m
    │   ├── plot_stress_T6.m
    │   ├── plot_stress_contour_Q4.m
    │   ├── plot_stress_contour_Q8.m
    │   ├── plot_stress_contour_T3.m
    │   ├── plot_stress_contour_T6.m
    │   └── print_displacement.m
    ├── Preprocess
    │   ├── create_2d_Q4_model.m
    │   ├── create_2d_Q8_model.m
    │   ├── create_2d_T3_model.m
    │   ├── create_2d_T6_model.m
    │   ├── read_fem_Q4_dat.m
    │   ├── read_fem_Q8_dat.m
    │   ├── read_fem_T3_dat.m
    │   ├── read_fem_T6_dat.m
    │   ├── read_fem_dat.asv
    │   └── setup_ID_LM.m
    ├── Solver
    │   └── solvedr.m
    └── Utils
        ├── assembly.m
        ├── gauss.m
        ├── gauss_2d.m
        ├── gauss_Tri3_2d.m
        └── plot_option.m


