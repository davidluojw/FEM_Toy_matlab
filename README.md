FEM_Toy_matlab
--------------

# Program Structure
```plaintext
FEM/
├── README.md
├── driver_2d_elasticity_linear_material.m
├── elasticity-mises.pdf
├── elasticity-sxx.pdf
├── input
│   ├── input_2d_Q4_16.dat
│   ├── input_2d_Q4_4.dat
│   ├── input_2d_Q4_testA.dat
│   ├── input_2d_Q4_testB.dat
│   ├── input_2d_Q4_testC.dat
│   ├── input_2d_T3_200.dat
│   ├── input_2d_T3_20000.dat
│   ├── input_2d_T3_50.dat
│   ├── input_2d_T3_500000.dat
│   └── input_2d_T3_8.dat
└── src
    ├── BoundaryCondition
    │   └── concentrate_traction.m
    ├── Element
    │   ├── Area_ShapeBasisN_2d.m
    │   ├── Area_ShapeBasisN_Grad_2d.m
    │   ├── Poly_ShapeBasisN_1d.m
    │   ├── Quad_ShapeBasisN_Grad_2d.m
    │   └── Quad_ShpaeBasisN_2d.m
    ├── Model
    │   ├── elasticity_elemQuad4_2d.m
    │   └── elasticity_elemTri3_2d.m
    ├── Postprocess
    │   ├── get_nodal_stress_Q4.m
    │   ├── get_nodal_stress_T3.m
    │   ├── get_qdpt_stress_Q4.m
    │   ├── get_qdpt_stress_T3.m
    │   ├── plot_displacement.m
    │   ├── plot_mesh.m
    │   ├── plot_stress_Q4.m
    │   ├── plot_stress_T3.m
    │   └── print_displacement.m
    ├── Preprocess
    │   ├── create_2d_Q4_model.m
    │   ├── create_2d_T3_model.m
    │   ├── read_fem_Q4_dat.m
    │   ├── read_fem_T3_dat.m
    │   ├── read_fem_dat.asv
    │   └── setup_ID_LM.m
    ├── Solver
    │   └── solvedr.m
    └── Utils
        ├── assembly.m
        ├── gauss.m
        ├── gauss_2d.m
        └── gauss_Tri3_2d.m


